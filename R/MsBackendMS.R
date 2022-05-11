#' @include hidden_aliases.R
NULL

#' @title MS data backend for ms files
#'
#' @aliases MsBackendMS-class
#'
#' @description
#'
#' The `MsBackendMS` class supports import of MS1/MS2 spectra data from
#' files in SIRIUS MS file format. `MsBackendMS` extends the
#' [MsBackendDataFrame()] backend directly and supports thus the
#' [applyProcessing()] function to make data manipulations persistent.
#'
#' New objects are created with the `MsBackendMS` function. The
#' `backendInitialize` method has to be subsequently called to
#' initialize the object and import MS data from (one or more) ms
#' files.
#'
#' The `MsBackendMS` backend provides an `export` method that allows to export
#' the data from the `Spectra` object (parameter `x`) to a file in ms format.
#' Parameters to this function are:
#'
#' - `x`: the `Spectra` object that should be exported.
#' - `file`: `character(1)` with the desired file name.
#' - `mapping`: named `character` providing the mapping between spectra
#'   variables and MS data fields. Defaults to
#'   `mapping = spectraVariableMapping(MsBackendMS())`.
#' - `allVariables`: `logical(1)` whether all spectra variables in `x` should be
#'   exported or only those defined with `mapping`.
#' - `exportName`: `logical(1)` whether a `NAME` field should always be exported
#'   even if not provided in `x`.
#' 
#' See the package vignette for details and examples.
#'
#' The `spectraVariableMapping` function allows to provide the mapping between
#' spectra variable names (i.e. the names that will be used for the spectra
#' variables in the [Spectra()] object) and the data field names of the
#' MS file. Parameter `format` allows to select pre-defined mappings. Currently
#' supported mapping flavors are:
#' 
#' - `format = "ms"`: default MS field names. Should work with SIRIUS 4.4.0
#'   files
#' 
#' @param object Instance of `MsBackendMS` class.
#'
#' @param file `character` with the (full) file name(s) of the ms file(s)
#'     from which MS/MS data should be imported or exported.
#'
#' @param format For `spectraVariableMapping`: `character(1)` specifying for
#'     which MS *flavour* the mapping should be returned. Currently supported
#'     are: `format = "ms"` (generic MS format)
#' 
#' @param mapping named `character` vector to rename MS fields to spectra
#'     variables (see [spectraVariableMapping()]). This allows to correctly
#'     import also custom fields or data from files with different MS
#'     *flavors*.
#' 
#' @param allVariables `logical(1)` whether all spectra variables in `x`
#'     should be exported or only those defined with `mapping`.
#'
#' @param exportName `logical(1)` whether a `NAME` field should always be
#'     exported even if not provided in `x`.
#' 
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup to import data in parallel. Defaults to `BPPARAM =
#'     bpparam()`. See [bpparam()] for more information.
#'
#' @param x For `export`: a [Spectra()] object that should be exported to the
#'     specified MS file.
#' 
#' @param ... Currently ignored.
#'
#' @return `MsBackendMS` and `backendInitialize` return an instance of a
#'     `MsBackendMS` class. `spectraVariableMapping` a named `character`
#'     vector with the mapping between spectra variables and MS data fields.
#' 
#' @author Steffen Neumann, Michael Witting, Laurent Gatto and Johannes Rainer
#'
#' @importClassesFrom Spectra MsBackendDataFrame
#'
#' @exportClass MsBackendMS
#'
#' @name MsBackendMS
#'
#' @examples
#'
#' ## Import spectra from a MS file 
#' f <- system.file("extdata", "example-sirius-440.ms",
#'     package = "MsBackendMS")
#' be <- backendInitialize(MsBackendMS(), f)
#' be
#'
#' be$msLevel
#' be$intensity
#' be$mz
#'
#' ## precursor m/z are however all missing
#' be$precursorMz
#'
#' ## Default spectra variable mapping
#' spectraVariableMapping(MsBackendMS())
#'
#'
#' ## Exporting Spectra objects in MS format.
#' 
#' sps <- Spectra(be)
#' export(MsBackendMS(), sps, file = stdout())
NULL

setClass("MsBackendMS",
         contains = "MsBackendDataFrame",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importMethodsFrom Spectra backendInitialize spectraData<- $<- $
#'
#' @importFrom BiocParallel bpparam
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @importFrom methods validObject
#'
#' @exportMethod backendInitialize
#'
#' @rdname MsBackendMS
setMethod("backendInitialize", signature = "MsBackendMS",
          function(object, file, 
                   mapping = spectraVariableMapping(object), ...,
                   BPPARAM = bpparam()) {
              if (missing(file) || !length(file))
                  stop("Parameter 'file' is mandatory for ", class(object))
              if (!is.character(file))
                  stop("Parameter 'file' is expected to be a character vector",
                       " with the file names from where data should be",
                       " imported")
              file <- normalizePath(file)
              if (any(!file.exists(file)))
                  stop("file(s) ",
                       paste(file[!file.exists(file)], collapse = ", "),
                       " not found")
              ## Import data and rbind.
              message("Start data import from ", length(file), " files ... ",
                      appendLF = FALSE)
              res <- bplapply(file, FUN = readMS, mapping = mapping,
                              BPPARAM = BPPARAM)
              message("done")
              res <- do.call(rbindFill, res)
              spectraData(object) <- res
              object$dataStorage <- "<memory>"
              validObject(object)
              object
          })

#' @rdname MsBackendMS
#'
#' @importFrom methods new
#'
#' @export MsBackendMS
MsBackendMS <- function() {
    new("MsBackendMS")
}

#' @importMethodsFrom Spectra spectraVariableMapping
#' 
#' @exportMethod spectraVariableMapping
#'
#' @rdname MsBackendMS
setMethod("spectraVariableMapping", "MsBackendMS",
          function(object, format = c("ms")) {
              switch(match.arg(format),
                     "ms" = c(
                         name = "NAME",
                         accession = "DB#",
                         formula = "FORMULA",
                         inchikey = "INCHIKEY",
                         adduct = "PRECURSORTYPE",
                         exactmass = "EXACTMASS",
                         rtime = "RETENTIONTIME",
                         precursorMz = "PRECURSORMZ",
                         adduct = "PRECURSORTYPE",
                         smiles = "SMILES",
                         inchi = "INCHI",
                         polarity = "IONMODE",
                         instrument = "INSTRUMENT"
                     )
                     )
          })

#' @importMethodsFrom Spectra export
#'
#' @exportMethod export
#'
#' @rdname MsBackendMS
setMethod("export", "MsBackendMS",
          function(object, x, file = tempfile(),
                   mapping = spectraVariableMapping(MsBackendMS()),
                   allVariables = TRUE, exportName = TRUE, ...) {
    if (missing(x))
        stop("Required parameter 'x' is missing. 'x' should be a 'Spectra' ",
             "object with the full spectra data.")
    if (!inherits(x, "Spectra"))
        stop("Parameter 'x' is supposed to be a 'Spectra' object with the full",
             " spectra data to be exported.")
    .export_ms(x = x, con = file, mapping = mapping,
                allVariables = allVariables, exportName = exportName)
})
