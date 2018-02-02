def get_gd(teff):
    if teff > 7250:
        # Radiative outer envelope
        return 1.0
    else:
        # Convective outer envelope
        return 0.32
