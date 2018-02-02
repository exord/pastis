class RVgaussianFitError(Exception):
    def __init__(self, mesg, contrast, rv0, sigma):
        self.mesg = mesg
        self.contrast = contrast
        self.rv0 = rv0
        self.sigma = sigma

    def __str__(self):
        return self.mesg


    
