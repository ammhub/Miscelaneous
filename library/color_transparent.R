# transparent color (stackoverflow)

    color.transparent <- function(color, a = 0.5){ #a[0,1]
        color_new <- col2rgb(color)
        apply(color_new, 2,
            function(x){
            rgb(red = x[1], green = x[2], blue = x[3], alpha = a*255, maxColorValue = 255)
            }
        )
    }
