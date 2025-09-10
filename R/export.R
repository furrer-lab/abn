#' Export abnFit object to structured format
#'
#' @param object An object of class abnFit
#' @param format The export format, currently only "json" is supported.
#' @param include_network Whether to include network structure
#' @param file Optional file path to save JSON. If NULL, returns JSON string
#' @param pretty Logical, whether to format JSON nicely
#' @param ... Additional export options
#' @return A JSON string or writes to file if 'file' is provided
export_abnFit <- function(object, format = "json", include_network = TRUE, file = NULL, pretty = TRUE, ...) {
  if (!inherits(object, "abnFit")) {
    stop("Input object must be of class 'abnFit'")
  }

  # Dispatch based on method
  if (object$method == "mle") {
    export_list <- export_abnFit_mle(object, format, include_network, ...)
  } else if (object$method == "bayes") {
    export_list <- export_abnFit_bayes(object, format, include_network, ...)
  } else {
    stop("Unsupported method in abnFit object. Supported methods are 'mle' and 'bayes'.")
  }

  # Convert to desired format
  if (format == "json") {
    return(export_to_json(export_list, format, file, pretty))
  } else {
    stop("Currently, only 'json' format is supported.")
  }
}

  json_str <- jsonlite::toJSON(json_list, pretty = pretty, auto_unbox = TRUE)

  if (!is.null(file)) {
    jsonlite::write_json(json_list, file, pretty = pretty, auto_unbox = TRUE)
    return(invisible(file))
  }

  return(json_str)
}
