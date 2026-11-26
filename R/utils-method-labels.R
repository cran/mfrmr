

.mfrm_model_code <- function(model) {
    model <- toupper(as.character(model[1L]))
    if (identical(model, "GMFRM")) model <- "GPCM"
    match.arg(model, c("RSM", "PCM", "GPCM"))
}
