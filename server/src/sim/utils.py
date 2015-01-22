def sanitize(arraywithnans):
        """ Sanitize input to remove NaNs. Warning, does not work on multidimensional data!! """
        from numpy import array, isnan
        arraywithnans = array(arraywithnans) # Make sure it's an array
        sanitized = arraywithnans[~isnan(arraywithnans)]
        if len(sanitized)==0:
                sanitized = 0
                print('                WARNING, no data entered for this parameter, assuming 0')

        return sanitized


def findinds(val1, val2=None, eps=1e-6):
    """
    Little function to find matches even if two things aren't eactly equal (eg. 
    due to floats vs. ints). If one argument, find nonzero values. With two arguments,
    check for equality using eps. Returns a tuple of arrays if val1 is multidimensional,
    else returns an array.
    
    Examples:
        findinds(rand(10)<0.5) # e.g. array([2, 4, 5, 9])
        findinds([2,3,6,3], 6) # e.g. array([2])
    
    Version: 2014nov27 by cliffk
    """
    from numpy import nonzero, array, ndim
    if val2==None: # Check for equality
        output = nonzero(val1) # If not, just check the truth condition
    else:
        output = nonzero(abs(array(val1)-val2)<eps) # If absolute difference between the two values is less than a certain amount
    if ndim(val1)==1: # Uni-dimensional
        output = output[0] # Return an array rather than a tuple of arrays if one-dimensional
    return output



def smoothinterp(newx, origx, origy, smoothness=10):
    """
    Smoothly interpolate over values and keep end points. Same format as numpy.interp.
    
    Example:
        from utils import smoothinterp
        origy = array([0,0.1,0.3,0.8,0.7,0.9,0.95,1])
        origx = linspace(0,1,len(origy))
        newx = linspace(0,1,5*len(origy))
        newy = smoothinterp(newx, origx, origy, smoothness=5)
        plot(newx,newy)
        hold(True)
        scatter(origx,origy)
    
    Version: 2014dec01 by cliffk
    """
    from numpy import interp, convolve, linspace, concatenate, ones, exp
    kernel = exp(-linspace(-2,2,2*smoothness+1)**2)
    kernel /= kernel.sum()
    newy = interp(newx, origx, origy) # Use interpolation
    newy = concatenate([newy[0]*ones(smoothness), newy, newy[-1]*ones(smoothness)])
    newy = convolve(newy, kernel, 'valid') # Smooth it out a bit
    return newy
    


def printarr(arr, arrformat='%0.2f  '):
    """ 
    Print a numpy array nicely.
    
    Example:
        from utils import printarr
        from numpy import random
        printarr(rand(3,7,4))
    
    Version: 2014dec01 by cliffk
    """
    from numpy import ndim
    if ndim(arr)==1:
        string = ''
        for i in range(len(arr)):
            string += arrformat % arr[i]
        print(string)
    elif ndim(arr)==2:
        for i in range(len(arr)):
            printarr(arr[i], arrformat)
    elif ndim(arr)==3:
        for i in range(len(arr)):
            print('='*len(arr[i][0])*len(arrformat % 1))
            for j in range(len(arr[i])):
                printarr(arr[i][j], arrformat)
    else:
        print(arr) # Give up
    return None