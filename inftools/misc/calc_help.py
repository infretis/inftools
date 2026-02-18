def estimate_error(arr, frac=0.5):
    """Estimates an error number based on last fraction
    of the running error.
    """
    import numpy as np

    if not 0 < frac <= 1:
        raise ValueError(f"frac must be between 0 and 1, got {frac}")
    arrl = int(len(arr)*frac)

    return np.average(arr[arrl:])

