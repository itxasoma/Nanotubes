# Magnetic Ordering in Chiral Nanotubes with Competing Interactions

The simulation consists of a FORTRAN Monte Carlo code that uses the Metropolis algorithm to flip one spin at a time, following the instructions of a three-dimensional (3D) classical Heisenberg Hamiltonian with exchange, dipolar and anisotropy interactions.

## Model

Our model is based on 3D Heisenberg spins $\vec{S}$ placed at the nodes of the nanotube lattice and interacting through the following Hamiltonian:

$$\hat{H} = E_{ex} + E_{dip} + E_{ani}$$

Where:

- **Exchange interaction**: 
$$E_{ex} = -\sum_{\langle i,j \rangle} J_{ij} \vec{S}_i \cdot \vec{S}_j$$

- **Dipolar interaction**: $E_{dip} = g \sum_{i<j} \frac{\vec{S}_i \cdot \vec{S}j - 3(\vec{S}i \cdot \hat{r}{ij})(\vec{S}j \cdot \hat{r}{ij})}{|\vec{r}{ij}|^3}$

- **Anisotropy**: 
$$E_{ani} = -K \sum_i (\vec{S}_i \cdot \vec{n}_i)^2$$

We define the dimensionless parameter $\gamma = g/J$ as the ratio between dipolar and exchange interactions.

## Cylindrical Components

To characterize the equilibrium configurations, we use the total magnetization components in cylindrical coordinates:

$$m_\rho = \frac{1}{N}\sum_{i=1}^N s_{\rho i}$$

$$m_\phi = \frac{1}{N}\sum_{i=1}^N |s_{\phi i}|$$

$$m_z = \frac{1}{N}\sum_{i=1}^N |s_{z i}|$$

where the cylindrical spin components are calculated as:

$$s_{\rho i} = \sin \theta_i \cos(\phi_i - \alpha_i)$$

$$s_{\phi i} = \sin \theta_i \sin(\phi_i - \alpha_i)$$  

$$s_{z i} = \cos \theta_i$$

## More details about the simulation code

### Methodology and computational details

For the tube construction, we simulate a 2D plane lattice with adjustable basis vectors and then roll it along the direction determined by the indices of the chiral vector to obtain the 3D nanotube. The tube geometry is defined by the pair $(N, N_z)$, where $N$ is the number of spins per layer and $N_z$ is the number of layers along the tube axis. The code also makes sure that the spin separation and tube radius are the same for both Zigzag and Angle configurations.

The code allocates arrays for spin vectors, neighbour lists, exchange, dipolar, and anisotropy fields. The constant parameters of the different interactions, $J$, $g$ and $K$ (uniaxial or radial) are specified in the input parameters.

A Monte Carlo loop is executed at each temperature point, first performing thermalization steps and then collecting averages of energy, magnetization, specific heat, and susceptibility. Temperature is updated following either a linear or annealing protocol. Observables are normalized and written to output files after averaging. Final spin configurations are saved in both cylindrical and planar projections for post-processing and visualization. 

Simulation results are stored in two types of files:
- Files named `Ener_N1_..._g_...out`, each containing 11 columns: Temperature, Total, Dipolar, Exchange, and Anisotropy Energies, $M_x$, $M_y$, $M_z$, $|\vec{M}|$, Specific Heat, and Susceptibility
- Files starting from `Conf_end`, which store spin vectors

`Writeaux` contains flags to control auxiliary output file generation. The file names `Pos`, `Config`, and `Nneigh` define the prefixes for geometry, spin, and neighbour data files, respectively.

The code is scalable and reproducible and can be applied to a wide range of scenarios. Simulations were run on CSUC's supercomputer *Pirineus III*.

### Simulation inputs

The simulation input allows customization of the nanotube structure and Monte Carlo parameters:

**Structure Parameters:**
- `Version`: selects either `Zigzag` or `Angle` construction method
- `Lattice`: defines the lattice type (`SC` for square or `General` for oblique geometries)
- `Alpha`: controls the angular relation between lattice vectors (e.g., $90°$ for square, $60°$ for triangular)
- `Natoms` and `Layers`: set the number of atoms per layer and number of layers respectively
- `n1`, `n2`: specify the chirality of the tube, determining its wrapping
- `nx`, `ny`: define the initial lattice dimensions before wrapping
- `Tlength`: scales the tube length along the axial direction

**Interaction Parameters:**
- `J_ex`: exchange interaction constant ($>0$ for ferromagnetic, $<0$ for antiferromagnetic coupling)
- `gdip`: sets the dipolar coupling strength constant ($g$), which scales with inverse cube of inter-spin distance
- `K_ani`: sets the anisotropy magnitude
- `V_ani`: selects the anisotropy axis (`Radial` or `Uniaxial`)

**Monte Carlo Parameters:**
- `MCtot` and `preMC`: define the number of Monte Carlo steps for measurement and thermalization
- `Protocol`: can be `Linear` or `Annealing`, controlling temperature evolution from `T_init` to `T_final`
- `dTemp`: step size or factor for temperature updates
- `New Orientation`: specifies trial move type (`Random`, `Ising`, or `Cone`)
- `Rmax`: defines the maximum cone aperture for trial orientations
- `SpinOrientation`: sets initial spin distribution (`Random`, `Alignedz`, `VortexZ`, etc.)

## Key Findings

The study reveals a rich variety of stable magnetic states including:
- **Uniform configurations**: Ferromagnetic (FM) and Antiferromagnetic (AFM) states
- **Vortex-like configurations**: Radial spin arrangements
- **Helical configurations**: Spiral magnetic ordering along the tube axis

The nature of these states strongly depends on:
- Tube aspect ratio and length
- Chirality (AA, AB, Zigzag types)
- The competition ratio $\gamma = g/J$
- Anisotropy effects
