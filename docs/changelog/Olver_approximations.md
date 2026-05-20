# Small-\(\chi\), shifted-\(q\), and Olver approximations for hyperspherical Bessel functions

This note summarizes approximations for the regular hyperspherical Bessel function
\(\phi_l(\chi)\), targeting errors measured as absolute error relative to the
peak amplitude of the function, rather than pointwise relative error.

Define

\[
S_K(\chi)=
\begin{cases}
\sin\chi, & K=1,\\
\chi, & K=0,\\
\sinh\chi, & K=-1.
\end{cases}
\]

The generalized Bessel equation being solved is

\[
\phi_l''(\chi)+2\frac{S_K'(\chi)}{S_K(\chi)}\phi_l'(\chi)
+\left[\nu^2-K-\frac{l(l+1)}{S_K^2(\chi)}\right]\phi_l(\chi)=0.
\]

Equivalently, defining the Schrödinger-like radial function

\[
\tilde u(\chi)=S_K(\chi)\phi_l(\chi),
\]

one obtains

\[
\tilde u''(\chi)+
\left[\nu^2-\frac{l(l+1)}{S_K^2(\chi)}\right]\tilde u(\chi)=0.
\]

For \(K=0\), the exact regular solution is

\[
\phi_l(\chi)=j_l(\nu\chi),
\qquad
\tilde u(\chi)=\chi j_l(\nu\chi).
\]

Use

\[
L=l(l+1),
\qquad
\ell=\sqrt{L}.
\]

For the full uniform Olver expression it is usually better to use the Langer
parameter

\[
\lambda=l+\frac12,
\]

whereas the local small-\(\chi\) and shifted-\(q\) approximations below use
\(\ell=\sqrt{l(l+1)}\).

For \(K=1\), all local criteria should use the folded value of \(\chi\) in the
standard symmetry interval, effectively \(0\leq\chi\leq\pi/2\).

---

## 1. Turning point

For a parameter \(\eta\), where \(\eta=\ell\) or \(\eta=\lambda\) depending on the
approximation, the turning point is defined by

\[
S_K(\chi_t)=\frac{\eta}{\nu}.
\]

When the approximation uses \(\ell\), this is

\[
\chi_t=
\begin{cases}
\asin(\ell/\nu), & K=1,\\
\ell/\nu, & K=0,\\
\asinh(\ell/\nu), & K=-1.
\end{cases}
\]

For the full Olver expression, replace \(\ell\) by \(\lambda=l+1/2\):

\[
\chi_t^{(\lambda)}=
\begin{cases}
\asin(\lambda/\nu), & K=1,\\
\lambda/\nu, & K=0,\\
\asinh(\lambda/\nu), & K=-1.
\end{cases}
\]

---

## 2. Full uniform Olver expression

Let

\[
\lambda=l+\frac12,
\qquad
\alpha=\frac{\nu}{\lambda},
\qquad
z_t=\frac{\lambda}{\nu}=\frac{1}{\alpha},
\qquad
S_K(\chi_t^{(\lambda)})=\frac{\lambda}{\nu}.
\]

The full leading-order Olver approximation writes the curved solution in terms
of a flat spherical Bessel function evaluated at an Olver coordinate \(z(\chi)\):

\[
\boxed{
\phi_l^{\rm Olver}(\chi)
\approx
\frac{z}{S_K(\chi)}
A(\chi)
 j_l(\nu z)
}
\]

where \(z=z(\chi)\) is determined by equality of the curved and flat actions,
measured from the turning point. In dimensionless form, define

\[
\mathcal I_K(\chi)=
\begin{cases}
\displaystyle
\int_\chi^{\chi_t}
\left[\frac{1}{S_K^2(s)}-\alpha^2\right]^{1/2}ds,
& \chi<\chi_t,\\[2ex]
\displaystyle
\int_{\chi_t}^{\chi}
\left[\alpha^2-\frac{1}{S_K^2(s)}\right]^{1/2}ds,
& \chi>\chi_t,
\end{cases}
\]

and

\[
\mathcal I_{\rm flat}(z)=
\begin{cases}
\displaystyle
\int_z^{z_t}
\left[\frac{1}{s^2}-\alpha^2\right]^{1/2}ds,
& z<z_t,\\[2ex]
\displaystyle
\int_{z_t}^{z}
\left[\alpha^2-\frac{1}{s^2}\right]^{1/2}ds,
& z>z_t.
\end{cases}
\]

The Olver coordinate is defined by

\[
\boxed{
\mathcal I_{\rm flat}(z)=\mathcal I_K(\chi)
}
\]

with \(z<z_t\) below the turning point and \(z>z_t\) above the turning point.

The amplitude is

\[
\boxed{
A(\chi)=
\left|
\frac{\alpha^2-z^{-2}}
     {\alpha^2-S_K^{-2}(\chi)}
\right|^{1/4}
}
\qquad (\chi\neq\chi_t).
\]

Equivalently,

\[
A(\chi)=
\left|
\frac{\nu^2-\lambda^2/z^2}
     {\nu^2-\lambda^2/S_K^2(\chi)}
\right|^{1/4}.
\]

At the turning point this has the finite limiting value

\[
A(\chi_t)=|S_K'(\chi_t)|^{-1/6},
\]

so

\[
A(\chi_t)=
\begin{cases}
\cos(\chi_t)^{-1/6}, & K=1,\\
1, & K=0,\\
\cosh(\chi_t)^{-1/6}, & K=-1.
\end{cases}
\]

### Analytic action formulae

Let

\[
y=\alpha S_K(\chi).
\]

The action is zero at \(y=1\). The branch \(y<1\) is below the turning point;
the branch \(y>1\) is above it. The following formulae use the continuous
branches of \(\log\), \(\acos\), and \(\atan2\), chosen so that
\(\mathcal I_K\geq0\) and \(\mathcal I_K=0\) at the turning point.

For the flat action, with \(y=\alpha z\),

\[
\mathcal I_{\rm flat}(y)=
\begin{cases}
\displaystyle
\log\left(\frac{1+\sqrt{1-y^2}}{y}\right)-\sqrt{1-y^2},
& y<1,\\[2ex]
\displaystyle
\sqrt{y^2-1}-\acos\left(\frac{1}{y}\right),
& y>1.
\end{cases}
\]

For \(K=-1\), i.e. \(S_K(\chi)=\sinh\chi\), define

\[
R_+=\sqrt{y^2-1},
\qquad
U=\sqrt{y^2+\alpha^2},
\]

for \(y>1\), and

\[
R_-=\sqrt{(1-y^2)(y^2+\alpha^2)}
\]

for \(y<1\). Then

\[
\mathcal I_{-1}(y)=
\begin{cases}
\displaystyle
\frac{\alpha}{2}
\log\left(
\frac{2y^2+\alpha^2-1+2R_+U}{1+\alpha^2}
\right)
-
\atan2(\alpha R_+,U),
& y>1,\\[3ex]
\displaystyle
\frac{\alpha}{2}
\atan2\left(-2R_-,2y^2+\alpha^2-1\right)
+
\frac12
\log\left(
\frac{2\alpha R_-+2\alpha^2+y^2(1-\alpha^2)}{y^2(1+\alpha^2)}
\right),
& y<1.
\end{cases}
\]

For \(K=1\), i.e. \(S_K(\chi)=\sin\chi\), define

\[
R_+=\sqrt{y^2-1},
\qquad
V=\sqrt{\alpha^2-y^2},
\]

for \(y>1\), and

\[
R_-=\sqrt{(1-y^2)(\alpha^2-y^2)},
\qquad
B=\frac{\alpha(1-y^2)}{R_-}
\]

for \(y<1\). Then

\[
\mathcal I_{1}(y)=
\begin{cases}
\displaystyle
\frac{\alpha}{2}
\atan2\left(2R_+V,\alpha^2+1-2y^2\right)
-
\atan2(\alpha R_+,V),
& y>1,\\[3ex]
\displaystyle
\frac12\log\left(\frac{1+B}{1-B}\right)
-
\frac{\alpha}{2}
\log\left(
\frac{\alpha^2-2y^2+1+2R_-}{\alpha^2-1}
\right),
& y<1.
\end{cases}
\]

The closed \(K=1\) formula assumes the folded interval where
\(0\leq\chi\leq\pi/2\), so \(0\leq y\leq\alpha\).

### Expected accuracy and range of validity

This is a uniform leading-order Olver approximation: it is designed to remain
smooth through the turning point. Unlike the local small-\(\chi\) approximations
below, it has no restriction of the form \(\chi\ll1\).

For \(K=0\), the expression reduces to the exact flat result. For \(K=\pm1\), the
accuracy is controlled by the neglected higher-order Liouville-Green terms, not
by a small-\(\chi\) truncation. It is expected to work best for moderate or large
\(l\) and \(\nu\), and should be treated cautiously for the very lowest
multipoles, where direct recursion is preferable. In practical use it is the
most robust of the approximations listed here, especially near or across the
turning point.

---

## 3. Local small-\(\chi\) Olver expansion

For the local expansion use

\[
\ell=\sqrt{l(l+1)},
\qquad
\alpha=\frac{\nu}{\ell},
\qquad
t=\alpha\chi=\frac{\nu\chi}{\ell},
\qquad
h=\frac{K}{\alpha^2}=K\frac{\ell^2}{\nu^2}.
\]

The differentiated Liouville-Green map is

\[
\left(\frac{dz}{d\chi}\right)^2
\left(\alpha^2-\frac{1}{z^2}\right)
=
\alpha^2-\frac{1}{S_K^2(\chi)}.
\]

Write

\[
z=\chi F(t,h),
\qquad
A=D(t,h)^{-1/2}.
\]

Through \(O(h^3)\),

\[
F=
1
-\frac{h}{6}
-\frac{h^2(4t^2+13)}{360}
-\frac{h^3(48t^4+148t^2+737)}{45360},
\]

and

\[
D=
1
-\frac{h}{6}
-\frac{h^2(12t^2+13)}{360}
-\frac{h^3(240t^4+444t^2+737)}{45360}.
\]

The local Olver approximation is therefore

\[
\boxed{
\phi_l^{\rm small\text{-}\chi}(\chi)
\approx
\frac{\chi F}{S_K(\chi)}
\frac{1}{\sqrt{D}}
 j_l(\nu\chi F)
}
\]

and should be rejected if

\[
D\leq0.
\]

### Error estimate

Let

\[
x=\frac{\nu\chi_{\max}}{\ell},
\qquad
r=|K|\frac{\ell^2}{\nu^2}.
\]

The first omitted terms in \(F\) and \(D\) give the estimates

\[
\epsilon_F=
 r^4
\frac{
576x^6+2256x^4+9432x^2+50801
}{5443200},
\]

and

\[
\epsilon_D=
 r^4
\frac{
4032x^6+11280x^4+28296x^2+50801
}{5443200}.
\]

A useful estimate of the error relative to the peak amplitude is

\[
\boxed{
\epsilon_{\rm small\text{-}\chi}
\approx
\nu\chi_{\max}\epsilon_F+\frac{1}{2}\epsilon_D.
}
\]

The first term estimates the Bessel phase error caused by the coordinate error
\(\delta z\). The second estimates the amplitude error from the missing terms in
\(D\).

A conservative target for \(10^{-4}\)-of-peak accuracy is

\[
\boxed{
\epsilon_{\rm small\text{-}\chi}<3\times10^{-5},
\qquad
D>0.
}
\]

A simpler conservative rule of thumb is

\[
\boxed{
\frac{\ell}{\nu}\lesssim0.55,
\qquad
\frac{\ell^2\chi_{\max}^7}{\nu}\lesssim0.3.
}
\]

A looser rule is

\[
\boxed{
\frac{\ell}{\nu}\lesssim0.6,
\qquad
\frac{\ell^2\chi_{\max}^7}{\nu}\lesssim1.
}
\]

The second condition follows from the large-\(x\) part of the first omitted term,
roughly

\[
\nu\chi_{\max}\epsilon_F
\sim
\frac{\ell^2\chi_{\max}^7}{9450\nu}.
\]

---

## 4. Amplitude-corrected shifted-\(q\) approximation

Expanding the curved centrifugal term gives

\[
\frac{1}{S_K^2(\chi)}
=
\frac{1}{\chi^2}
+\frac{K}{3}
+\frac{K^2\chi^2}{15}
+\frac{2K^3\chi^4}{189}
+\cdots.
\]

Keeping only the constant curvature correction gives

\[
q^2=\nu^2-\frac{K\ell^2}{3}.
\]

The amplitude-corrected shifted approximation is

\[
\boxed{
\phi_l^{q}(\chi)
\approx
\left(\frac{q}{\nu}\right)^{1/2}
\frac{\chi}{S_K(\chi)}
 j_l(q\chi)
}
\]

or equivalently

\[
\boxed{
\phi_l^{q}(\chi)
\approx
\left(1-\frac{K\ell^2}{3\nu^2}\right)^{1/4}
\frac{\chi}{S_K(\chi)}
 j_l(q\chi).
}
\]

The amplitude prefactor is important because

\[
\left(\frac{q}{\nu}\right)^{1/2}
=
\left(1-\frac{K\ell^2}{3\nu^2}\right)^{1/4}
=
1-\frac{K\ell^2}{12\nu^2}+O\left(\frac{\ell^4}{\nu^4}\right),
\]

matching the leading amplitude correction in the local Olver expansion.

The approximation requires

\[
\boxed{q^2>0.}
\]

For \(K=-1\), this is automatic. For \(K=1\), it requires

\[
\nu^2>\frac{\ell^2}{3}.
\]

### Refined shifted-\(q\) error proxy

The older separated estimates

\[
\frac{\ell^2\chi_{\max}^3}{90\nu},
\qquad
\frac{1}{180}
\left(\frac{\ell}{\nu}\right)^4
\left[2+\left(\frac{\nu\chi_{\max}}{\ell}\right)^2\right]
\]

are useful for scaling, but they are not well calibrated across the range of
\(\chi_{\max}\). In particular, simply replacing a \(3\times10^{-5}\) cut by
\(3\times10^{-4}\) tends to be too restrictive at larger \(\chi_{\max}\), while
being too permissive around \(\chi_{\max}\simeq0.2\)--\(0.3\).

A better practical proxy is obtained by comparing the amplitude-corrected
shifted formula directly with the local Olver small-\(\chi\) expansion.
Let

\[
\alpha=\frac{\nu}{\ell},
\qquad
t_{\max}=\alpha\chi_{\max},
\qquad
|h|=\alpha^{-2}.
\]

The shifted argument corresponds to

\[
\frac{q}{\nu}
=
\left(1-\frac{h}{3}\right)^{1/2}.
\]

Comparing this with the Olver coordinate \(z=\chi F\) gives the leading argument
mismatch

\[
\Delta(\nu z)
\approx
\frac{\ell\,t_{\max}(t_{\max}^2+2)}{90\alpha^4}.
\]

The previous cubic phase estimate only kept the \(t_{\max}^3\) part; the
additional \(2t_{\max}\) term is important at smaller \(\chi_{\max}\).

For the amplitude mismatch, compare

\[
\frac{F}{\sqrt D}
\]

with

\[
\left(\frac{q}{\nu}\right)^{1/2}
=
\left(1-\frac{h}{3}\right)^{1/4}.
\]

The leading residual amplitude mismatch is

\[
\Delta A
\approx
\frac{|t_{\max}^2-2|}{180\alpha^4}.
\]

This improves on the envelope \((2+t_{\max}^2)/(180\alpha^4)\), which is too
pessimistic at larger \(t_{\max}\) and misses the cancellation near
\(t_{\max}\simeq\sqrt2\).

Define

\[
\boxed{
\epsilon_{\rm shift}
=
 c_{\rm ph}
\frac{\ell\,t_{\max}(t_{\max}^2+2)}{90\alpha^4}
+
\frac{|t_{\max}^2-2|}{180\alpha^4}.
}
\]

The coefficient \(c_{\rm ph}\) converts raw argument mismatch into a peak-amplitude
error proxy. A useful empirical value is

\[
\boxed{c_{\rm ph}\simeq0.5,}
\]

while \(c_{\rm ph}=1\) is a more conservative analytic choice.

For a target peak-amplitude error of about \(3\times10^{-4}\), a reasonable
calibrated cut is

\[
\boxed{
q^2>0,
\qquad
\epsilon_{\rm shift}\lesssim 7\times10^{-4}.
}
\]

Empirically, for a slightly looser target of about \(4\times10^{-4}\), the simple
cut

\[
\boxed{
q^2>0,
\qquad
\epsilon_{\rm shift}<10^{-3}
}
\]

appears sufficient in practice.

### Turning-point check

If the interval reaches the turning point,

\[
\chi_{\max}\gtrsim\chi_t,
\]

where

\[
\chi_t=
\begin{cases}
\asin(\ell/\nu), & K=1,\\
\ell/\nu, & K=0,\\
\asinh(\ell/\nu), & K=-1,
\end{cases}
\]

then also check that the shifted Bessel turning point agrees with the curved
one:

\[
\boxed{
|q\chi_t-\ell|\lesssim\hbox{a few}\times10^{-4}.
}
\]

For small \(\ell/\nu\),

\[
q\chi_t-\ell
\simeq
\frac{\ell^5}{30\nu^4}.
\]

In practice this turning-point condition is usually implied by the refined
\(\epsilon_{\rm shift}\) cut whenever \(\chi_{\max}\gtrsim\chi_t\).

---

## 5. Recommended hierarchy

For \(K=0\), use the exact result

\[
\phi_l(\chi)=j_l(\nu\chi).
\]

For \(K\neq0\):

1. Use the amplitude-corrected shifted-\(q\) approximation only in the small domain
   where

   \[
   q^2>0,
   \qquad
   \epsilon_{\rm shift}\lesssim7\times10^{-4}
   \]

   for a target error of about \(3\times10^{-4}\). For a practical target around
   \(4\times10^{-4}\), the empirical condition

   \[
   q^2>0,
   \qquad
   \epsilon_{\rm shift}<10^{-3}
   \]

   is often sufficient.

2. Otherwise use the local small-\(\chi\) Olver expansion if

   \[
   D>0,
   \qquad
   \epsilon_{\rm small\text{-}\chi}<3\times10^{-5},
   \]

   or, more simply, if

   \[
   \frac{\ell}{\nu}\lesssim0.55,
   \qquad
   \frac{\ell^2\chi_{\max}^7}{\nu}\lesssim0.3.
   \]

3. Otherwise use the full uniform Olver approximation. It is the best of these
   approximations near the turning point and outside the small-\(\chi\) domain.
