from numpy import dtype, float64

# TODO: cleanup and add doc

flow_type = dtype([("rho",   float64),
                   ("u",     float64),
                   ("v",     float64),
                   ("w",     float64),
                   ("E",     float64),
                   ("p",     float64),
                   ("H",     float64),
                   ("T",     float64),
                   ("mu",    float64),
                   ("k",     float64),
                   ("c",     float64),
                   ("cp",    float64),
                   ("gamma", float64)])
