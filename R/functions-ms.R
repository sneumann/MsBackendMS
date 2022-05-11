#' @title Reading MS files
#'
#' The `readMS` function imports the data from a file in *.ms format used by e.g. SIRIUS software reading
#' all specified fields and returning the data as a [DataFrame()].
#'
#' @param f `character(1)` with the path to an MS file.
#'
#' @param msLevel `numeric(1)` with the MS level. Default is 2. This value will
#'     be reported as the spectra's MS level **unless** the source MS file
#'     defines the MS level.
#'
#' @param mapping named `character` vector to rename MS fields to spectra
#'     variables (see [spectraVariableMapping()]). This allows to correctly
#'     import also custom fields or data from files with different MS
#'     *flavors*.
#'
#' @param ... Additional parameters, currently ignored.
#'
#' @return
#'
#' A `DataFrame` with each row containing the data from one spectrum
#' in the MS file. m/z and intensity values are available in columns `"mz"`
#' and `"intensity"` in a list representation.
#'
#' @export
#'
#' @importFrom Spectra coreSpectraVariables
#' 
#' @importFrom S4Vectors DataFrame
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom MsCoreUtils rbindFill
#'
#' @importFrom stats ave
#' 
#' @importFrom methods as
#'
#' @author Laurent Gatto, Steffen Neumann, Johannes Rainer
#'
#' @examples
#'
#' fls <- dir(system.file("extdata", package = "MsBackendMS"),
#'     full.names = TRUE, pattern = "ms$")[1L]
#'
#' readMsp(fls)
readMS <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMS()), ...) {
    if (length(f) != 1L)
        stop("Please provide a single ms file.")
    
    ms <- scan(file = f, what = "",
                sep = "\n", quote = "",
                allowEscapes = FALSE,
                quiet = TRUE)
    
    ## Find individual records
    begin <- grep("^>compound", ms, ignore.case = FALSE)
    end <- c(begin[-1] -1L, length(ms))

    sp <- mapply(begin, end, FUN = function(a, b) {
         .extract_ms_spectrum(ms[a:b])
    }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    res <- DataFrame(rbindFill(sp))

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res$mz <- NumericList(res$mz, compress = FALSE)
    res$intensity <- NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if (!any(colnames(res) == msLevel))
        res$msLevel <- as.integer(msLevel)
    res
}

#' @param ms `character()` of lines defining a spectrum in ms
#'     format.
#'
#' @param mapping spectra variable mapping that allows renaming data fields.
#' 
#' @author Laurent Gatto, Johannes Rainer
#' 
#' @importFrom stats setNames
#'
#' @noRd
.extract_ms_spectrum <- function(ms) {
    ## grep metadata
    desc.idx <- grep("^>", ms)
    spec.idx <- grep("^>ms", ms)
    
    ## Remove duplicates from description, implies removing spectra beginnings:
    desc.idx <- desc.idx[ave(c(desc.idx, spec.idx), c(desc.idx, spec.idx), FUN = length) == 1]

    ms1 <- .extract_ms_peaks(ms[(spec.idx[1]+1):(spec.idx[2]-1)])
    ms2 <- .extract_ms_peaks(ms[(spec.idx[2]+1):length(ms)])
    
    if (!length(ms1)) 
        ms1 <- matrix(numeric(), ncol = 2L)
    else if (is.unsorted(ms1[, 1L]))
        ms1 <- ms1[order(ms1[, 1L]), ]

    if (!length(ms2)) 
        ms2 <- matrix(numeric(), ncol = 2L)
    else if (is.unsorted(ms2[, 1L]))
        ms2 <- ms2[order(ms2[, 1L]), ]
    
    ## Handle tag/Value
    desc <- ms[desc.idx]
    r <- regexpr("^>(?<tag>[^ ]+) (?<value>.+)", desc, perl=TRUE)

    desc <- .parse.tagvalue(desc, r)
    

    # ## map fields to spectra variables
    # idx <- match(names(desc), mapping)
    # not_na <- !is.na(idx)
    # if (any(not_na))
    #     names(desc)[not_na] <- names(mapping)[idx][not_na]
    # 
    # ## Handle eventually duplicated names -> list
    # if (anyDuplicated(names(desc))) {
    #     res <- split(unname(desc), names(desc))
    #     dups <- lengths(res) > 1L
    #     dup_res <- res[dups]
    #     res <- as.data.frame(res[!dups])
    #     for (name in names(dup_res))
    #         res <- do.call("$<-", list(res, name, unname(dup_res[name])))
    # } else res <- as.data.frame(as.list(desc))
    
    res <- as.data.frame(as.list(desc))
    
    ## Ensure correct data type
    ## polarity
    if (any(have <- colnames(res) == "polarity"))
        res[, have] <- .process_polarity(res[, have])
    if (any(have <- colnames(res) == "msLevel"))
        res[, have] <- .process_mslevel(res[, have])

    res$mz <- list(ms2[, 1L])
    res$intensity <- list(ms2[, 2L])
    res
}

.extract_ms_peaks <- function(ms) {
    res <- do.call(rbind, strsplit(sub("^(\\t|[[:space:]]+)", "", ms),
                                  "[[:space:]]+"))
    mode(res) <- "double"
    res
}

## from regexr manpage
.parse.tagvalue <- function(res, result) {
    m <- do.call(rbind, lapply(seq_along(res), function(i) {
        if(result[i] == -1) return("")
        st <- attr(result, "capture.start")[i, ]
        substring(res[i], st, st + attr(result, "capture.length")[i, ] - 1)
    }))
    colnames(m) <- attr(result, "capture.names")
    setNames(m[,2], m[,1])
}


    
.process_polarity <- function(x, input = TRUE) {
    if (input) {
        if (grepl("^(p|\\+)", x, ignore.case = TRUE))
            return(1L)
        if (grepl("^(n|-)", x, ignore.case = TRUE))
            return(0L)
        -1L
    } else {
        xnew <- rep(NA_character_, length(x))
        xnew[x == 1L] <- "Positive"
        xnew[x == 0L] <- "Negative"
        xnew
    }
}

#' @param x value to be formatted
#'
#' @param input `logical(1)` whether the data is imported or exported.
#'
#' @noRd
.process_mslevel <- function(x, input = TRUE) {
    if (input)
        as.integer(sub("ms", "", x, ignore.case = TRUE))
    else paste0("MS", x)
}

#' @description
#'
#' Function to export a `Spectra` object in MS format to `con`.
#'
#' @param x `Spectra`
#'
#' @param con output file.
#'
#' @param mapping named `character` vector that maps from `spectraVariables`
#'    (i.e. `names(mapping)`) to the variable name that should be used in the
#'    MS file.
#'
#' @param allVariables `logical(1)` whether all spectra variables in `x` should
#'    be exported or only those that are listed in `mapping`. Note that if
#'    `exportName = TRUE` a field *NAME* will be exported regardless of
#'    `mapping`.
#'
#' @param exportName `logical(1)` whether a NAME field will always be exported
#'    even if no such spectra variable is available in `x`.
#' 
#' @author Michael Witting, Johannes Rainer
#'
#' @importMethodsFrom Spectra spectraVariables spectraNames
#'
#' @importMethodsFrom Spectra peaksData spectraData
#'
#' @noRd
#'
#' @examples
#'
#' spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
#' spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
#' spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
#'
#' sps <- Spectra(spd)
#'
#' .export_ms(sps)
#'
#' ## Handling of variables with multiple entries
#' sps$synonym <- list(c("a", "b"), "d", c("e", "f", "g"))
#' .export_ms(sps)
#'
.export_ms <- function(x, con = stdout(),
                        mapping = spectraVariableMapping(MsBackendMS()),
                        allVariables = TRUE, exportName = TRUE) {
    spv <- spectraVariables(x)
    spv <- spv[!(spv %in% c("dataOrigin", "dataStorage"))]
    if (!allVariables)
        spv <- spv[spv %in% names(mapping)]
    spd <- spectraData(x, spv)
    ## Process any known required data conversions
    if (any(spv == "msLevel"))
        spd$msLevel <- .process_mslevel(spd$msLevel, input = FALSE)
    if (any(spv == "polarity"))
        spd$polarity <- .process_polarity(spd$polarity, input = FALSE)
    idx <- match(colnames(spd), names(mapping))
    colnames(spd)[!is.na(idx)] <- mapping[idx[!is.na(idx)]]
    ## Force variable NAME:
    if (!any(tolower(colnames(spd)) == "name") && exportName)
        spd$NAME <- seq_len(nrow(spd))
    idx <- which(tolower(colnames(spd)) == "name")
    if (length(idx)) {
        idx <- idx[1L]
        spd <- spd[, c(idx, seq_len(ncol(spd))[-idx])]
    }

    ## Determine which columns contain list-like data (i.e. multiple entries).
    mult <- colnames(spd)[!vapply(spd, function(z) is.vector(z) & !is.list(z),
                                  logical(1))]
    for (m in mult) {
        spd[, m] <- vapply(
            spd[, m], function(z) paste0(z, collapse = paste0("\n", m, ": ")),
            character(1))
    }

    tmp <- lapply(colnames(spd), function(z) {
        paste0(z, ": ", spd[, z], "\n")
    })

    pks <- vapply(peaksData(x), function(z)
        paste0("Num Peaks: ", nrow(z), "\n",
               paste0(paste0(z[, 1], " ", z[, 2], "\n"), collapse = ""),
               collapse = ""),
        character(1))
    tmp <- do.call(cbind, c(tmp, list(pks)))
    tmp[grep(": NA\n", tmp, fixed = TRUE)] <- ""
    writeLines(apply(tmp, 1, paste0, collapse = ""), con = con)
}
