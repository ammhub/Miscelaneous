# color rectangle area in a plot - graphics (stackoverflow)

    color.rectangle <- function(color){ 
        params <- par()$usr 
        rect(params[1], params[3], params[2], params[4], col = color) 
    } 
