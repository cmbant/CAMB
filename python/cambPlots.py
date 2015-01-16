# Handy routines from plotting CAMB/CAMB_sources outputs
# See covCompare.py for sample script to plot outs

from pylab import *
import numpy as np

def_colors = ['b', 'r', 'm', 'c', 'g', 'k', 'y', '0.25', '0.75']

def saveCls(f, ls, vecs):
    np.savetxt(f, hstack([ls.reshape(-1, 1) ] + [x.reshape(-1, 1) for x in vecs]), fmt='%5i ' + '%12.7e ' * len(vecs))

def cl_plot(ax, L, cl, x=None, y=None, grid=False, lmax=None, **axargs):
    if lmax is None: lmax = int(L[-1])
    lmax_map = lmax
    if all(cl == 0): return
    if x == 'sqrt':
        ls = np.sqrt(L)
        lmax_map = sqrt(lmax)
    else: ls = L[:]
    handle = ax.plot(ls, cl, **axargs)
    if x == 'sqrt':
#            ticks = [2, 10, 40, 100, 200, 400, 700]
        ticks = [2, 40, 200, 700]
        if (lmax >= 1500): ticks += range(1500, lmax, 1000)
        ax.set_xticks(sqrt(ticks))
        ax.set_xticklabels(ticks)
        print ls
    elif x == 'log': ax.set_xscale('log')
    if y == 'log': ax.set_yscale('log')
    ax.set_xlim(ls[0], lmax_map)
    if grid: ax.grid(True)
    return handle


class ClResult(object):

    def __init__(self, lensedcl=None, lens_potential=None, cov_cls=None, unlensedcl=None, scaling=1, skiprows=0):
        if isinstance(lensedcl, basestring): lensedcl = np.loadtxt(lensedcl, skiprows=skiprows)
        if lensedcl is not None:
            self.set_lensedcl(lensedcl, scaling)
        if isinstance(lens_potential, basestring): lens_potential = np.loadtxt(lens_potential, skiprows=skiprows)
        if lens_potential is not None:
            self.set_lensingcl(lens_potential, scaling)
        if isinstance(unlensedcl, basestring): unlensedcl = np.loadtxt(unlensedcl, skiprows=skiprows)
        if unlensedcl is not None:
            self.set_unlensedcl(unlensedcl, scaling)

        if cov_cls is not None:
            self.set_cov_cls(cov_cls, scaling)

    def set_cov_cls(self, cov_cls, scaling=1, tag_names=None):
        arr = np.loadtxt(cov_cls)
        n = np.int(np.sqrt(arr.shape[1] - 1))
        self.l = np.rint(arr[:, 0])
        for i in range(n):
            arr[:, 1 + 2 + i * n] *= sqrt(self.l * (self.l + 1))  # deflection angle rather than lensing potential
            arr[:, 1 + i + n * 2] *= sqrt(self.l * (self.l + 1))  # deflection angle rather than lensing potential

        self.lmax = self.l[-1]
        self.cl_array = [[arr[:, x + y * n + 1] for x in range(n)] for y in range(n)]
        self.cl_windows = [[arr[:, x + y * n + 1] for x in range(3, n)] for y in range(3, n)]
        if tag_names is None:
            self.cl_array_tags = ['T', 'E', '$\\phi$'] + [('$w_' + str(i + 1) + '$') for i in range(n - 3)]
        else: self.cl_array_tags = tag_names
        self.cl_window_tags = [self.cl_array_tags[i] for i in range(3, n)]

    def set_lensedcl(self, lensedcl, scaling=1):
        self.l = np.rint(lensedcl[:, 0])
        self.TT = lensedcl[:, 1] * scaling
        self.EE = lensedcl[:, 2] * scaling
        self.BB = lensedcl[:, 3] * scaling
        self.TE = lensedcl[:, 4] * scaling
        self.cls = lensedcl[:, 1:5] * scaling
        self.cl_array = [[self.TT, self.TE, None], [self.TE, self.EE, None], [None, None, self.BB]]
        self.cl_array_tags = ['T', 'E', 'B']
        self.lmax = self.l[-1]

    def set_unlensedcl(self, cl, scaling=1):
        # all internal CL are in TT, EE, BB, TE order
        self.unlens_cls = np.zeros((cl.shape[0], 4))
        self.unlens_cls[:, 0:2] = cl[:, 1:3] * scaling
        self.unlens_cls[:, 3] = cl[:, 4] * scaling
        self.unlensl = np.rint(cl[:, 0])

    def set_lensingcl(self, cl, scaling=1):
        # lens_potential_output_file
        self.unlens_cls = cl[:, 1:5] * scaling
        self.lensl = np.rint(cl[:, 0])
        self.unlensl = self.lensl
        self.phiphi = cl[:, 5]
        self.phiT = cl[:, 6] * sqrt(scaling)
        self.phiE = cl[:, 7] * sqrt(scaling)
        self.phicls = cl[:, 5:8]
        self.phicls[:, 1:2] *= sqrt(scaling)
        self.lmax_phi = self.lensl[-1]


    def compare_array_acc(self, C, lmax=None, lmin=2, fails=None,
                      tols=None, array_tag='cl_array', tags=None):
        if hasattr(self, array_tag) and hasattr(C, array_tag):
            c1 = getattr(self, array_tag)
            c2 = getattr(C, array_tag)
            imin = int(lmin - self.l[0])
            if lmax is None: lmax = min(self.lmax, C.lmax) - 50
            imax = int(lmax - self.l[0])
            if tols is None: tols = [1e-3, 3e-3, 3e-3]
            for i, _ in enumerate(self.cl_array):
                for i2 in range(i):
                    if c1[i][i2] is not None:
                        diff = c1[i][i2][imin:imax] - c2[i][i2][imin:imax]
                        if i != i2:
                            norm = np.sqrt(c1[i][i][imin:imax] * c1[i2][i2][imin:imax] + c1[i][i2][imin:imax] ** 2)
                        else: norm = c1[i][imin:imax]
                        diff /= norm
                        if np.any(abs(diff) > sqrt(tols[i] * tols[i2])):
                            rms = np.sqrt(np.average(diff ** 2))
                            if tags is None: tags = self.cl_array_tags
                            label = tags[i] + tags[i2]
                            if fails is None: fails = dict()
                            plot(self.l[imin:imax], diff)
                            fails[label] = rms
        return fails

    def compare_type_tol(self, C, tag='TT', lmax=None, lmin=2, fails=None, norm=None, tol=1e-3):
        if hasattr(self, tag) and hasattr(C, tag):
            c1 = getattr(self, tag)
            c2 = getattr(C, tag)
            imin = lmin - self.l[0]
            if lmax is None: lmax = min(self.lmax, C.lmax)
            imax = lmax - self.l[0]
            if norm is None: norm = c1
            diff = (c1[imin:imax] - c2[imin:imax]) / norm[imin:imax]
            if np.any(abs(diff) > tol):
                rms = np.sqrt(np.average(diff ** 2))
                if fails is None: fails = dict()
                plot(self.l[imin:imax], diff)
                fails[tag] = rms

        return fails

    def compare_lensing(self, C, **pars):
        return self.compare_type(C, 'phiphi', **pars)



def checkClResult(result, lensedCl=True):
    if isinstance(result, basestring):
        if lensedCl: return ClResult(result)
        else: return ClResult(lens_potential=result)
    elif isinstance(result, np.ndarray): return ClResult(result)
    if not isinstance(result, ClResult): raise Exception('expecting a ClResult object, numpy array or string for filename')
    return result

def plot_array(clResults, lmax=None, lmin=2, indices=None, array_tag='cl_array', tags=None, diff=False,
                   label=True, half=False, figsize=(12, 8), tag_separator='', colors=['k', 'r', 'b', 'g', 'm'], **axargs):
        if isinstance(clResults, ClResult): clResults = [clResults]
        if isinstance(clResults, basestring): clResults = [clResults]
        for i, res in enumerate(clResults):
            if isinstance(res, basestring): clResults[i] = ClResult(cov_cls=res)

        R1 = clResults[0]
        imin = int(lmin - R1.l[0])
        if lmax is None: lmax = min([R.lmax for R in clResults]) - 50
        if indices is None: indices = range(len(getattr(R1, array_tag)))
        f, plots = subplots(len(indices), len(indices), figsize=figsize, sharex='col')
        plots = np.reshape(plots, (len(indices), len(indices)))  # in case only 1
        imax = int(lmax - R1.l[0])
        for ir, result in enumerate(clResults):
            c1 = getattr(result, array_tag)
            if diff and ir == 0:
                compare = c1
            else:
                for ix, i in enumerate(indices):
                    for ix2, i2 in enumerate(indices):
                        if not half or i2 <= i:
                            if c1[i][i2] is not None:
                                data = c1[i][i2][imin:imax]
                                if diff:
                                    comp = compare[i][i2][imin:imax]
                                    delta = data - comp
                                    if i != i2:
                                        norm = np.sqrt(compare[i][i][imin:imax] * compare[i2][i2][imin:imax] + comp ** 2)
                                    else: norm = comp
                                    data = delta / norm

                                if tags is None: tags = R1.cl_array_tags
                                cl_plot(plots[ix, ix2], R1.l[imin:imax], data, color=colors[ir], **axargs)
                                if label:
                                    text(0.9, 0.9, tags[i] + tag_separator + tags[i2], horizontalalignment='right',
                                         verticalalignment='top', transform=plots[ix, ix2].transAxes)
                        else: plots[ix, ix2].axis('off')
        tight_layout()
        return f, plots


def plot_compare(clResults, colors=def_colors, ncl=1, x='sqrt', x_lens='log', y=None, y_lens='log', compare=None, diff_fraction=True, **args):
    if isinstance(clResults, ClResult): clResults = [clResults]
    if isinstance(clResults, basestring): clResults = [clResults]
    for i, res in enumerate(clResults):
        clResults[i] = checkClResult(res)
    rcParams['figure.figsize'] = max(8, 5 * ncl), 8
    f, ax = plt.subplots(2, ncl, sharex='col')
#        f.set_size_inches(max(8,5*ncl),8)
    ax = np.reshape(ax, (2, ncl))
    if compare is not None: cl_base = checkClResult(compare, lensedCl=True)
    for i, cl in enumerate(clResults):
        if i == 0 and compare is None: cl_base = cl
        for icl in range(ncl):
            if icl <= 3:
                if (all(cl.cls[:, icl] == 0)): continue
                cl_plot(ax[0, icl], cl.l, cl.cls[:, icl], color=colors[i], x=x, y=y, **args)
                if i == 0 and compare is not None:
                    cl_plot(ax[0, icl], cl_base.l, cl_base.cls[:, icl], color='k', ls='--', x=x, **args)
                if i > 0 or compare is not None:
                    nmax = min(cl.l.shape[0], cl_base.l.shape[0])
                    if 'lmax' in args: nmax = min(nmax, args['lmax'] - 1)
                    diff = cl.cls[:nmax, icl] - cl_base.cls[:nmax, icl]
                    if icl < 3:
                        if diff_fraction: diff /= cl_base.cls[:nmax, icl]
                    elif icl == 3:
                        if diff_fraction: diff /= sqrt(cl_base.TT[:nmax] * cl_base.EE[:nmax])
                    cl_plot(ax[1, icl], cl.l[:nmax], diff, color=colors[i], y=None, x=x, **args)
            else:  # lensing
                cl_plot(ax[0, icl], cl.lensl, cl.phiphi, color=colors[i], x=x_lens, y=y_lens, **args)
                if i > 0 or compare is not None:
                    nmax = min(cl.lensl.shape[0], cl_base.lensl.shape[0])
                    diff = (cl.phiphi[:nmax] - cl_base.phiphi[:nmax]) / cl_base.phiphi[:nmax]
                    cl_plot(ax[1, icl], cl.lensl[:nmax], np.abs(diff), color=colors[i], x=x_lens, y=y_lens, **args)
    f.tight_layout()



