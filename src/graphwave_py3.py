import math
import networkx as nx
import numpy as np
import scipy as sc
import os

TAUS = [1, 10, 25, 50]
ORDER = 30
PROC = 'approximate'
ETA_MAX = 0.95
ETA_MIN = 0.80
NB_FILTERS = 2


def laplacian(a):
    n_nodes, _ = a.shape
    posinv = np.vectorize(lambda x: float(1.0) / np.sqrt(x) if x > 1e-10 else 0.0)
    d = sc.sparse.diags(np.array(posinv(a.sum(0))).reshape([-1, ]), 0)
    lap = sc.sparse.eye(n_nodes) - d.dot(a.dot(d))
    return lap
def charac_function(time_points, temp):
    temp2 = temp.T.tolil()
    d = temp2.data
    n_timepnts = len(time_points)
    n_nodes = temp.shape[1]
    final_sig = np.zeros((2 * n_timepnts, n_nodes))
    zeros_vec = np.array([1.0 / n_nodes*(n_nodes - len(d[i])) for i in range(n_nodes)])
    for i in range(n_nodes):
        final_sig[::2, i] = zeros_vec[i] +\
                            1.0 / n_nodes *\
                            np.cos(np.einsum("i,j-> ij",
                                             time_points,
                                             np.array(d[i]))).sum(1)
    for it_t, t in enumerate(time_points):
        final_sig[it_t * 2 + 1, :] = 1.0 / n_nodes * ((t*temp).sin().sum(0))

    return final_sig


def charac_function_multiscale(heat, time_points):
    final_sig = []
    for i in heat.keys():
        final_sig.append(charac_function(time_points, heat[i]))
    return np.vstack(final_sig).T

def compute_cheb_coeff_basis(scale, order):
    xx = np.array([np.cos((2 * i - 1) * 1.0 / (2 * order) * math.pi)
                   for i in range(1, order + 1)])
    basis = [np.ones((1, order)), np.array(xx)]
    for k in range(order + 1-2):
        basis.append(2* np.multiply(xx, basis[-1]) - basis[-2])
    basis = np.vstack(basis)
    f = np.exp(-scale * (xx + 1))
    products = np.einsum("j,ij->ij", f, basis)
    coeffs = 2.0 / order * products.sum(1)
    coeffs[0] = coeffs[0] / 2
    return list(coeffs)


def heat_diffusion_ind(graph, taus=TAUS, order = ORDER, proc = PROC):
    '''
    This method computes the heat diffusion waves for each of the nodes
    INPUT:
    -----------------------
    graph    :    Graph (etworkx)
    taus     :    list of scales for the wavelets. The higher the tau,
                  the better the spread of the heat over the graph
    order    :    order of the polynomial approximation
    proc     :    which procedure to compute the signatures (approximate == that
                  is, with Chebychev approx -- or exact)
    OUTPUT:
    -----------------------
    heat     :     tensor of length  len(tau) x n_nodes x n_nodes
                   where heat[tau,:,u] is the wavelet for node u
                   at scale tau
    taus     :     the associated scales
    '''
    # Compute Laplacian
    a = nx.adjacency_matrix(graph)
    n_nodes, _ = a.shape
    thres = np.vectorize(lambda x : x if x > 1e-4 * 1.0 / n_nodes else 0)
    lap = laplacian(a)
    n_filters = len(taus)
    if proc == 'exact':
        ### Compute the exact signature
        lamb, U = np.linalg.eigh(lap.todense())
        heat = {}
        for i in range(n_filters):
             heat[i] = U.dot(np.diagflat(np.exp(- taus[i] * lamb).flatten())).dot(U.T)
    else:
        heat = {i: sc.sparse.csc_matrix((n_nodes, n_nodes)) for i in range(n_filters) }
        monome = {0: sc.sparse.eye(n_nodes), 1: lap - sc.sparse.eye(n_nodes)}
        for k in range(2, order + 1):
             monome[k] = 2 * (lap - sc.sparse.eye(n_nodes)).dot(monome[k-1]) - monome[k - 2]
        for i in range(n_filters):
            coeffs = compute_cheb_coeff_basis(taus[i], order)
            heat[i] = np.sum([ coeffs[k] * monome[k] for k in range(0, order + 1)]) # sc.sum to np.sum for new version
            temp = thres(heat[i].A) # cleans up the small coefficients
            heat[i] = sc.sparse.csc_matrix(temp)
    return heat, taus


def graphwave_alg(graph, time_pnts, taus= 'auto',
              verbose=False, approximate_lambda=True,
              order=ORDER, proc=PROC, nb_filters=NB_FILTERS,
              **kwargs):
    ''' wrapper function for computing the structural signatures using GraphWave
    INPUT
    --------------------------------------------------------------------------------------
    graph             :   nx Graph
    time_pnts         :   time points at which to evaluate the characteristic function
    taus              :   list of scales that we are interested in. Alternatively,
                          'auto' for the automatic version of GraphWave
    verbose           :   the algorithm prints some of the hidden parameters
                          as it goes along
    approximate_lambda:   (boolean) should the range oflambda be approximated or
                          computed?
    proc              :   which procedure to compute the signatures (approximate == that
                          is, with Chebychev approx -- or exact)
    nb_filters        :   nuber of taus that we require if  taus=='auto'
    OUTPUT
    --------------------------------------------------------------------------------------
    chi               :  embedding of the function in Euclidean space
    heat_print        :  returns the actual embeddings of the nodes
    taus              :  returns the list of scales used.
    '''
    if taus == 'auto':
        if approximate_lambda is not True:
            a = nx.adjacency_matrix(graph)
            lap = laplacian(a)
            try:
                l1 = np.sort(sc.sparse.linalg.eigsh(lap, 2,  which='SM',return_eigenvectors=False))[1]
            except:
                l1 = np.sort(sc.sparse.linalg.eigsh(lap, 5,  which='SM',return_eigenvectors=False))[1]
        else:
            l1 = 1.0 / graph.number_of_nodes()
        smax = -np.log(ETA_MIN) * np.sqrt( 0.5 / l1)
        smin = -np.log(ETA_MAX) * np.sqrt( 0.5 / l1)
        taus = np.linspace(smin, smax, nb_filters)
    heat_print, _ = heat_diffusion_ind(graph, list(taus), order=order, proc = proc)
    chi = charac_function_multiscale(heat_print, time_pnts)
    return chi, heat_print, taus

