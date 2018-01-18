import numpy as np
import matplotlib.pyplot as plt

def sinc_interp(d,x,xi):
    """ Perform sinc interpolation on a sampled signal.
    
    inputs
        d  -- uniformly sampled data points, spaced by 1
        x  -- uniform sample locations
        xi -- locations to evaluation for the sinc interpolation
    
    outputs
        di -- sinc interpolated values at locations xi
    """

    # Get d, x looking the way we expect with singleton dims
    x = np.expand_dims(x,axis=0)
    d = np.expand_dims(d,axis=0)
    
    X = 1 # always spaced by 1
    s = np.sinc((xi - np.transpose(x))/X)
    di = np.dot(d,s)

    return(np.squeeze(di))

def sinc_resample(dn,xn,xu):
    """ Perform sinc interpolation on a nonuniformly sampled signal.
    
    inputs
        dn -- non-uniformly sampled data points
        xn -- non-uniform sample locations
        xu -- uniform sample points, spaced by 1

    outputs
        du -- uniformly sampled data
    """

    # Get d,x looking the way we expect with singleton dims
    xu = np.expand_dims(xu,axis=0)
    
    X = 1 # always spaced by 1
    E = np.sinc((xn - np.transpose(xu))/X)
    du = np.linalg.lstsq(np.transpose(E),dn.flatten())[0]

    return(du)

if __name__ == '__main__':

    ## Introduction
    d = np.zeros((1,10))
    d = np.append(d,range(10,0,-1))
    d = np.append(d,range(11))
    d = np.append(d,np.zeros((1,10)))
    x = range(-20,21)

    fig1 = plt.figure()
    ax = fig1.add_subplot(211)
    plt.stem(x,d)
    plt.xlabel('x')
    plt.ylabel('s(x)')
    plt.draw()


    ## (1) Sinc Interpolation
    # Take the test signal, upsample by factor of 10
    xi = np.linspace(-20,20,401,endpoint=True)
    di = sinc_interp(d,x,xi);

    ax = fig1.add_subplot(212)
    plt.plot(xi,di,'k')
    plt.stem(x,d)
    plt.xlabel('x')
    plt.ylabel('s(x)')
    plt.title('Original Data and sinc interpolation')
    plt.draw()


    ## (2) Bunched Samples
    # Generate the bunched sample data by subsampling the sinc interpolated
    # data
    db = di[2::20]
    db = np.append(db,di[7::20])
    xb = xi[2::20]
    xb = np.append(xb,xi[7::20])

    fig2 = plt.figure()
    ax = fig2.add_subplot(211)
    plt.plot(xi,di,label='Original Sinc Interpolation')
    plt.stem(xb,db,label='Bunched Samples')
    plt.xlabel('x')
    plt.ylabel('s(x)')
    plt.title('Bunched Samples')
    ax.legend()
    plt.draw()

    # Now see if we can recover the original samples from the bunched samples
    du = sinc_resample(db,xb,x)
    ax = fig2.add_subplot(212)
    plt.plot(x,du,label='Recovered')
    plt.plot(x,d,'--',label='Original')
    plt.stem(xb,db,label='Bunched Samples')
    plt.xlabel('x')
    plt.ylabel('s(x)')
    plt.title('Bunched Samples and sinc resampling')
    ax.legend()
    plt.draw()

    ## (3) Random Samples
    # Choose random samples from the sinc interpolated signal
    ndx = np.random.choice(range(1,len(di)),size=41,replace=False)
    dr = di[ndx]
    xr = xi[ndx]

    fig3 = plt.figure()
    ax = fig3.add_subplot(211)
    plt.stem(xr,dr)
    plt.title('Random Samples')
    plt.draw()

    dur = sinc_resample(dr,xr,x)
    ax = fig3.add_subplot(212)
    plt.plot(x,dur,label ='Recovered')
    plt.plot(x,d,'--',label='Original')
    plt.stem(xr,dr,label='Random Samples')
    ax.legend()
    plt.title('Random Samples and sinc resampling')
    plt.draw()
    
    # Show all the plots
    plt.show()
