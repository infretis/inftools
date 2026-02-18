def estimate_error(arr, frac=0.5):
    """Estimates an error number based on last fraction
    of the running error.

    if frac is 0.2, average uses the last 20% of arr data.
    """
    import numpy as np

    # ensure frac is between 0 and 1
    if not 0 < frac <= 1:
        raise ValueError(f"frac must be between 0 and 1, got {frac}")

    # find location of where to start
    arrl = int(len(arr)*(1-frac))

    return np.average(arr[arrl:])

