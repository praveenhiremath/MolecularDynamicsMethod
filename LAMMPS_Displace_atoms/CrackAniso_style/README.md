# Implementation of Anisotropic displacement fields in LAMMPS package (A popular package for Molecular Dynamics method) 

The stress intensity-controlled loading was applied through utilization of the mode I displacement field, which is  proportional to the stress intensity factor, $K_I$, and can be obtained from anisotropic linear elastic fracture mechanics as \cite{sih1965cracks,sih1968fracture,cheung1994molecular}:
\begin{multline}
u_x=\frac{K_I \sqrt{2r}}{\sqrt{\pi}} \Re \bigg\{\Big[\frac{1}{(\mu_1-\mu_2)}\Big]\Big[\mu_1 p_2(\cos\theta+\mu_2 \sin\theta)^{1/2}\\
    -\mu_2 p_1(\cos\theta+\mu_1 \sin\theta)^{1/2}\Big] \bigg\}  ,  
\label{eq:uxx}    
\end{multline} 

\begin{multline}
u_y=\frac{K_I \sqrt{2r}}{\sqrt{\pi}} \Re\bigg\{\Big[\frac{1}{(\mu_1-\mu_2)}\Big]\Big[\mu_1 q_2(\cos\theta+\mu_2 \sin\theta)^{1/2}\\
    -\mu_2 q_1(\cos\theta+\mu_1 \sin\theta)^{1/2}\Big] \bigg\} , 
\label{eq:uyy}    
\end{multline}

\noindent where $\Re$ represents the real part operator, and
\begin{multline}
p_1=s_{11}\mu_1^2 +s_{12}-s_{16}\mu_1   \mbox{,  } p_2=s_{11}\mu_2^2 +s_{12}-s_{16}\mu_2,  \\
q_1=\frac{s_{12}\mu_1^2 +s_{22}-s_{26}\mu_1}{\mu_1} \mbox{,  } q_2=\frac{s_{12}\mu_2^2 +s_{22}-s_{26}\mu_2}{\mu_2}.
\end{multline}

\noindent The complex valued $\mu_1$ and $\mu_2$ were obtained by solving the characteristic equation
\begin{equation}
    s_{11}\mu_j ^4-2s_{16}\mu_j ^3+(2s_{12}+s_{66})\mu_j ^2-2s_{26}\mu_j+s_{22}=0,
\end{equation}
where the compliance constants $s_{ij}$, which were evaluated at $T=0$ K, were rotated to match the orientation of interest through suitable rotation operations.

