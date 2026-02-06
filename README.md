# Newton's Second Law For Rotation

We know that by Newton's second law, forces and momenta are tied
together by the following equation: $$\begin{equation}
    \boldsymbol{F}= \frac{d\boldsymbol{p}}{dt}
\end{equation}$$ When forces are balanced, this leads to a conservation
law (i.e. of linear momentum) and translational equilibrium. A similar
thing happens with a so-called term angular momentum (this time,
however, resulting in rotational equilibrium), which is defined below:
$$\begin{equation}
    \boldsymbol{L}= \boldsymbol{r}\times\boldsymbol{p}
\end{equation}$$ Here, $\boldsymbol{p}$ is the linear momentum of some
particle and $\boldsymbol{r}$ is its position vector. Hence, obviously,
this term's value is dependent on the choice of coordinate system.

Now, we define our "rotational force". The torque due to force
$\boldsymbol{F}$ which acts on a particle with position vector
$\boldsymbol{r}$ is defined as: $$\begin{equation}
    \boldsymbol{\tau}= \boldsymbol{r}\times\boldsymbol{F}
\end{equation}$$ Yet again, our choice of origin affects the value of
the torque. Now, for our equation that relates these two terms:
$$\begin{equation}
    \frac{d\boldsymbol{L}}{dt} = \frac{d}{dt}\left(\boldsymbol{r}\times\boldsymbol{p}\right) = \left(\frac{d\boldsymbol{r}}{dt}\times\boldsymbol{p}\right) + \left(\boldsymbol{r}\times\frac{d\boldsymbol{p}}{dt}\right) = \left(\frac{d\boldsymbol{r}}{dt}\times m\frac{d\boldsymbol{r}}{dt}\right) + \left(\boldsymbol{r}\times\boldsymbol{F}\right)
\end{equation}$$ $$\begin{equation}
    \boldsymbol{\tau}= \frac{d\boldsymbol{L}}{dt}
\end{equation}$$

# Fixed Axis Rotation

Now, say you have a body rotating about some axis. WLOG, let's call this
the z-axis. Further, assume that the origin will be along the z-axis.
Consider particle $j$ on this body and its angular momentum:
$$\begin{equation}
    \boldsymbol{L}(j) = \boldsymbol{r}_j\times m_j\boldsymbol{v}_j
\end{equation}$$ The path traced by this particle when revolving around
a fixed axis will always be circular, regardless of the object's shape.
So, by circular motion, we know that: $$\begin{equation}
    \boldsymbol{v}_j = \boldsymbol{\omega}\times\boldsymbol{r}_j
\end{equation}$$ The angular velocity points along the z-axis.
Additionally, the position vector will have a z-axis component and a
component on the xy-plane. The former would vanish. The latter would be
the perpendicular distance from the fixed axis to the particle,
$\rho_j$. Hence, the particle's velocity lies in the xy-plane and has
magnitude $$\begin{equation}
    v_j = \omega \rho_j
\end{equation}$$ We want to calculate the angular momentum about the
fixed-axis, z. We can find the dot product between $\boldsymbol{L}(j)$
and $\hat{z}$: $$\begin{equation}
    L_z(j) = \hat{z}\cdot (\boldsymbol{r}_j\times m_j\boldsymbol{v}_j)
\end{equation}$$ We can, again, break down $\boldsymbol{r}_j$ into a
component along the z-axis and on the xy-plane. The cross-product of the
former with $m_j\boldsymbol{v}_j$ would not align with the z-axis.
Hence, it would vanish. Thus, all that would remain would be the xy
component of $\boldsymbol{r}_j$. This would be perpendicular to
$\boldsymbol{v}_j$ and will have magnitude $\rho_j$. Hence:
$$\begin{equation}
    L_z(j) = m_j\rho_j^2\omega
\end{equation}$$ Thus, the total angular momentum of the body about the
z-axis with an origin aligning with this axis is: $$\begin{equation}
    L_z = \sum_j m_j\rho_j^2\omega = \omega\sum_j m_j\rho_j^2
\end{equation}$$ The summation term is called the moment of inertia of
the body. It is defined by: $$\begin{equation}
    I_{body} = \sum_jm_j\rho_j^2
\end{equation}$$ Then: $$\begin{equation}
    L_z = I_{body}\omega
\end{equation}$$ Naturally, the moment of inertia would be solved using
using an integral. A more useful definition would be: $$\begin{equation}
    I = \int\rho^2dm
\end{equation}$$ This can be calculated using some basic calculus for a
variety of different objects. The dynamic consequence of this is
straightforward: $$\begin{equation}
    \boldsymbol{\tau}= \frac{d}{dt}(I\omega) = I\alpha
\end{equation}$$ Here, $\alpha$ is the angular acceleration of the body.
This is the $F = ma$ for rotational bodies.

# Tensor Of Inertia

When rotating about a fixed axis, we found that the angular momentum
about that axis is equal to $I\omega$. When we eliminate the constraint
of a fixed axis, we find deeper complexities.

Let us begin by the most general expression for the angular momentum:
$$\begin{equation}
    \boldsymbol{L}= \sum_{j=1}^N(\boldsymbol{r}_j\times m_j\dot{\boldsymbol{r}}_j)
\end{equation}$$ Let's define the center of mass of a body as the
weighted average position (with the weight being associated with mass):
$$\begin{equation}
    \boldsymbol{R}= \frac{\sum_{j=1}^N m_j\boldsymbol{r}_j}{M}
\end{equation}$$ This position will turn out to be special and will
simplify the calculations. This is because, as we shall see, all that
matters to describe the rotation of a body is the angular momentum about
its center of mass. So, let's try writing (3.1) in terms of
$\boldsymbol{R}$: $$\begin{equation}
    \boldsymbol{r}_j = \boldsymbol{R}+ \boldsymbol{r}_j'
\end{equation}$$ $$\begin{equation}
    \dot{\boldsymbol{r}}_j = \dot{\boldsymbol{R}} + \dot{\boldsymbol{r}}_j'
\end{equation}$$ Here, $\dot{\boldsymbol{r}}_j'$ is the position of the
particle relative to the COM. Next: $$\begin{equation}
    \boldsymbol{L}= \sum_{j=1}^N\left(\left(\boldsymbol{R}+ \boldsymbol{r}_j'\right)\times m_j\left(\dot{\boldsymbol{R}} + \dot{\boldsymbol{r}}_j'\right)\right) = \sum_{j=1}^Nm_j\left(\left(\boldsymbol{R}+ \boldsymbol{r}_j'\right)\times\left(\dot{\boldsymbol{R}} + \dot{\boldsymbol{r}}_j'\right)\right)
\end{equation}$$ $$\begin{equation}
    = \sum_{j=1}^Nm_j\left(\boldsymbol{R}\times\dot{\boldsymbol{R}} + \boldsymbol{R}\times\dot{\boldsymbol{r}}_j' + \boldsymbol{r}_j'\times\dot{\boldsymbol{R}} + \boldsymbol{r}_j'\times \dot{\boldsymbol{r}}_j'\right)
\end{equation}$$ $$\begin{equation}
    =\boldsymbol{R}\times\sum m_j\dot{\boldsymbol{R}} + \boldsymbol{R}\times\sum m_j\dot{\boldsymbol{r}}_j' + \sum m_j\boldsymbol{r}_j'\times\dot{\boldsymbol{R}} + \sum m_j\left(\boldsymbol{r}_j'\times\dot{\boldsymbol{r}}_j'\right)
\end{equation}$$ It can be easily shown that
$\sum m_j\boldsymbol{r}_j' = \boldsymbol{0}$. This would also mean that
its derivative, $\sum m_j\dot{\boldsymbol{r}}_j' = 0$. Hence:
$$\begin{equation}
    \boldsymbol{L}= \boldsymbol{R}\times\sum m_j\dot{\boldsymbol{R}} + \sum m_j\left(\boldsymbol{r}_j'\times\dot{\boldsymbol{r}}_j'\right)
\end{equation}$$ $$\begin{equation}
    \boldsymbol{L}= \left(\boldsymbol{R}\times\dot{\boldsymbol{R}}\sum m_j\right) + \left(\sum \boldsymbol{r}_j'\times m_j\dot{\boldsymbol{r}}_j'\right) = \boldsymbol{R}\times M\boldsymbol{V}+ \sum \boldsymbol{r}_j'\times m_j\dot{\boldsymbol{r}}_j'
\end{equation}$$ One half of the angular momentum comes from the
translational motion of the COM. The other half comes from rotation
about the COM. Similarly, we can come up an expression for the torque on
the object: $$\begin{equation}
    \boldsymbol{\tau}= \sum\boldsymbol{r}_j\times\boldsymbol{F}_j = \sum\left(\boldsymbol{r}_j' + \boldsymbol{R}\right)\times\boldsymbol{F}_j = \sum\boldsymbol{r}_j'\times\boldsymbol{F}_j + \sum\boldsymbol{R}\times\boldsymbol{F}_j
\end{equation}$$ $$\begin{equation}
=\sum\boldsymbol{r}_j'\times\boldsymbol{F}_j + \boldsymbol{R}\times\sum\boldsymbol{F}_j = \boldsymbol{R}\times\boldsymbol{F}+ \sum\boldsymbol{r}_j'\times\boldsymbol{F}_j
\end{equation}$$ Relating (3.11) with (3.9) using (1.5) we get:
$$\begin{equation}
    \boldsymbol{R}\times\boldsymbol{F}+ \sum\boldsymbol{r}_j'\times\boldsymbol{F}_j = \frac{d}{dt}\left(\boldsymbol{R}\times M\boldsymbol{V}\right) + \frac{d}{dt}\left(\boldsymbol{r}_j'\times m_j\dot{\boldsymbol{r}}_j'\right)
\end{equation}$$ $$\begin{equation}
    \boldsymbol{R}\times\boldsymbol{F}+ \sum\boldsymbol{r}_j'\times\boldsymbol{F}_j = \boldsymbol{V}\times M\boldsymbol{V}+ \boldsymbol{R}\times M\boldsymbol{A}+ \frac{d}{dt}\left(\boldsymbol{r}_j'\times m_j\dot{\boldsymbol{r}}_j'\right)
\end{equation}$$ $$\begin{equation}
\sum\boldsymbol{r}_j'\times\boldsymbol{F}_j = \frac{d}{dt}\left(\boldsymbol{r}_j'\times m_j\dot{\boldsymbol{r}}_j'\right)
\end{equation}$$ This proves that center of mass's motion is irrelevant
and that: $$\begin{equation}
    \boldsymbol{\tau}_{COM} = \frac{d\boldsymbol{L}_{COM}}{dt}
\end{equation}$$ The angular momentum about a body's COM was:
$$\begin{equation}
    \boldsymbol{L}_{COM} = \sum\boldsymbol{r}_j'\times m_j\dot{\boldsymbol{r}}_j' = \sum\boldsymbol{r}_j'\times m_j\left(\boldsymbol{\omega}\times\boldsymbol{r}_j'\right)
\end{equation}$$ These cross products can be computed and each component
of the angular momentum can be found. The results are summarized below:
$$\begin{equation}
    L_x = \sum m_j\left(y_j^2 + z_j^2\right)\omega_x - \sum m_j x_j y_j \omega_y - \sum m_j x_j z_j\omega_z
\end{equation}$$ $$\begin{equation}
    L_y = \sum m_j\left(z_j^2 + x_j^2\right)\omega_y - \sum m_j y_j z_j \omega_z - \sum m_j y_j x_j\omega_x
\end{equation}$$ $$\begin{equation}
    L_z = \sum m_j\left(x_j^2 + y_j^2\right)\omega_z - \sum m_j z_j x_j \omega_x - \sum m_j z_j y_j\omega_y
\end{equation}$$ You may already see where this is going. The angular
momentum of a rigid body can be concisely written in terms of its
angular velocity, $\boldsymbol{\omega}$: $$\begin{equation}
    \boldsymbol{L}= \boldsymbol{I}\boldsymbol{\omega}
\end{equation}$$

# Kinetic Energy Of A Rigid Body

Later on, we will apply the Lagrangian to solve the rigid body dynamics
of the T-handle. This will, first, require the computation of kinetic
energy. $$\begin{equation}
    T = \frac{\sum m_j v_j^2}{2}
\end{equation}$$ Using equation (3.4), we get: $$\begin{equation}
    T = \frac{\sum m_j\left|\boldsymbol{V}+ \boldsymbol{v}_j'\right|^2}{2} = \frac{\sum m_j\left(\boldsymbol{V}+ \boldsymbol{v}_j'\right)\cdot\left(\boldsymbol{V}+ \boldsymbol{v}_j'\right)}{2}
\end{equation}$$ $$\begin{equation}
    = \frac{\sum m_j\left(V^2 + 2\boldsymbol{V}\cdot\boldsymbol{v}_j' + \left(v_j'\right)^2\right)}{2} = \frac{\sum m_jV^2 + 2\boldsymbol{V}\cdot\sum m_j\boldsymbol{v}_j' + \sum m_j\left(v_j'\right)^2}{2}
\end{equation}$$ The middle term cancels out by the fact that
$\sum m_jv_j'$ is zero. This yields the formula: $$\begin{equation}
    T=\frac{1}{2}MV^2 + \frac{1}{2}\sum m_j\left(v_j'\right)^2
\end{equation}$$ Next, let's process the rotational kinetic energy term:
$$\begin{equation}
    T_{rot} = \frac{1}{2}\sum m_j\left(v_j'\right)^2 = \frac{1}{2}\sum m_j\left(\boldsymbol{\omega}\times\boldsymbol{r}_j'\right)\cdot\left(\boldsymbol{\omega}\times\boldsymbol{r}_j'\right)
\end{equation}$$ A linear algebra trick can be used to convert this to:
$$\begin{equation}
    T_{rot} = \frac{1}{2}\sum m_j\boldsymbol{\omega}\cdot\left(\boldsymbol{r}_j'\times(\boldsymbol{\omega}\times\boldsymbol{r}_j')\right) = \frac{\boldsymbol{\omega}}{2}\cdot\sum m_j\boldsymbol{r}_j'\times(\boldsymbol{\omega}\times\boldsymbol{r}_j') = \frac{\boldsymbol{\omega}}{2}\cdot\sum\boldsymbol{r}_j'\times m_j(\boldsymbol{\omega}\times\boldsymbol{r}_j')
\end{equation}$$ $$\begin{equation}
    T_{rot} = \frac{\boldsymbol{\omega}\cdot\boldsymbol{L}_{COM}}{2} = \frac{\boldsymbol{\omega}\cdot\boldsymbol{I}\boldsymbol{\omega}}{2}
\end{equation}$$ For any body, it is possible to find a set of
coordinate axes such that any cross terms in the inertia matrix (such as
$I_{xy}$) equal to zero. These are known as principal axes.[^1] Then,
the inertia tensor becomes diagonal: $$\begin{equation}
    T_{rot} = \frac{1}{2}
    \begin{bmatrix}
        w_x\\
        w_y\\
        w_z
    \end{bmatrix}
    \cdot
    \begin{bmatrix}
    I_{xx}w_x\\
    I_{yy}w_y\\
    I_{zz}w_z
    \end{bmatrix} = \frac{1}{2}\left(I_{xx}w_x^2 + I_{yy}w_y^2 + I_{zz}w_z^2\right)
\end{equation}$$ $$\begin{equation}
    T = \frac{1}{2}MV^2 + \frac{L_x^2}{2I_{xx}} + \frac{L_y^2}{2I_{yy}} + \frac{L_z^2}{2I_{zz}}
\end{equation}$$

# Euler Rigid Body Equations

We have been using an inertial frame of reference up to this point.
However, recall that to diagonalize our inertia tensor, our coordinate
axes must align with the principal axes. This would mean that the
coordinate axes must rotate with the body, making it a non-inertial
frame. To apply (1.5), we must take the time derivative of the angular
momentum in an inertial frame. But, to simplify the inertia tensor, we
must be in a non-inertial frame (i.e. body frame). This is the
fundamental problem.

To resolve this issue, we must find the relation between the derivative
of angular momentum in both frames.

Let the angular momentum of a body about its COM in terms of the body
frame's basis vectors be defined by the following vector:
$$\begin{equation}
    \boldsymbol{L}_{COM} = L_x\hat{b}_x + L_y\hat{b}_y + L_z\hat{b}_z
\end{equation}$$ Taking the derivative of this in the inertial frame
yields: $$\begin{equation}
    \frac{d\boldsymbol{L}_{COM}}{dt} = \frac{d}{dt}\left(\sum L_j\hat{b}_j\right) = \left(\sum \dot{L}_j\hat{b}_j\right) + \left(\sum L_j\frac{d\hat{b}_j}{dt}\right)
\end{equation}$$ The LHS is simply the derivative in the inertial frame.
The RHS has two components. One is the derivative in the body frame. The
other shows how the body frame changes with respect to an inertial
frame. $$\begin{equation}
    \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{i} = \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{b} + \left(\sum L_j\frac{d\hat{b}_j}{dt}\right) = \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{b} + \left(\sum L_j\boldsymbol{\omega}\times\hat{b}_j\right)
\end{equation}$$ $$\begin{equation}
    = \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{b} + \boldsymbol{\omega}\times\sum L_j\hat{b}_j = \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{b} + \boldsymbol{\omega}\times\boldsymbol{L}_{COM}
\end{equation}$$ The body frame derivative can be evaluated:
$$\begin{equation}
    \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{b} = \frac{d}{dt}\Bigr|_{b} \left(L_x\hat{b}_x + L_y\hat{b}_y + L_z\hat{b}_z\right)
\end{equation}$$ In the body's frame of reference, the basis vectors are
stationary. Thus, they can be pulled out, giving the following result:
$$\begin{equation}
    \frac{d\boldsymbol{L}_{COM}}{dt}\Bigr|_{b} = \begin{bmatrix}
        \dot{L}_x\\
        \dot{L}_y\\
        \dot{L}_z
    \end{bmatrix}
\end{equation}$$ The term
$\boldsymbol{\omega}\times\boldsymbol{L}_{COM}$ can be evaluated by:
$$\begin{equation}
    \boldsymbol{\omega}\times\boldsymbol{L}_{COM} = \begin{bmatrix}
        \hat{i}&\hat{j}&\hat{k}\\
        w_x&w_y&w_z\\
        L_x&L_y&L_z
    \end{bmatrix} = \begin{bmatrix}
        w_yL_z - L_yw_z\\
        w_zL_x - L_zw_x\\
        w_xL_y - L_xw_y
    \end{bmatrix}
\end{equation}$$ Hence, RHS of (5.4) becomes: $$\begin{equation}
    \begin{bmatrix}
        w_yL_z - L_yw_z + \dot{L}_x\\
        w_zL_x - L_zw_x + \dot{L}_y\\
        w_xL_y - L_xw_y + \dot{L}_z
    \end{bmatrix} = \begin{bmatrix}
        w_yI_{zz}w_z - I_{yy}w_yw_z + I_{xx}\dot{w}_x\\
        w_zI_{xx}w_x - I_{zz}w_zw_x + I_{yy}\dot{w}_y\\
        w_xI_{yy}w_y - I_{xx}w_xw_y + I_{zz}\dot{w}_z
    \end{bmatrix}
\end{equation}$$ $$\begin{equation}
    \begin{bmatrix}
        w_yI_{zz}w_z - I_{yy}w_yw_z + I_{xx}\dot{w}_x\\
        w_zI_{xx}w_x - I_{zz}w_zw_x + I_{yy}\dot{w}_y\\
        w_xI_{yy}w_y - I_{xx}w_xw_y + I_{zz}\dot{w}_z
    \end{bmatrix} = \begin{bmatrix}
        \tau_x\\
        \tau_y\\
        \tau_z
    \end{bmatrix}
\end{equation}$$ The results are the Euler Rigid Body Equations. We
typically replace $I_{xx}$ with $I_x$ for all axes when principal axes
are involved because there are no cross terms. Thus:

$$\begin{equation}
    \dot{w}_x = \frac{\tau_x + (I_y - I_z)w_yw_z}{I_x}
\end{equation}$$ $$\begin{equation}
    \dot{w}_y = \frac{\tau_y + (I_z - I_x)w_zw_x}{I_y}
\end{equation}$$ $$\begin{equation}
    \dot{w}_z = \frac{\tau_z + (I_x - I_y)w_xw_y}{I_z}
\end{equation}$$ Note that these equations are in the rigid body frame.
dd

# Dzhanibekov Simulation

Note that

[^1]: We are using a body frame in this case.
