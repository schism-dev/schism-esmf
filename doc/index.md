# SCHISM-ESMF 

The Earth System Modeling Framework (ESMF) and National Unified Operational Prediction Capability (NUOPC) caps for SCHISM.

## Three steps

The cap is meant to be used in the context of ESMF-based coupling frameworks, which come in two flavors, either based on ESMF or on its NUOPC extension. We also provide some standalone example applications which you can use to base your own implementation on

1. Install the prerequisites following [Installation](installation.html)
2. Build the cap's ESMF and NUOPC libraries following [Installation](installation.html)
3. Include the SCHISM-ESMF cap in your application's build system by linking to `-libschism_esmf`, see below

## Coupling frameworks

### NUOPC

The National Unified Operational Prediction Capability (NUOPC) is a consortium of Navy, NOAA, and Air Force modelers and their research partners. It aims to advance the weather prediction modeling systems used by meteorologists, mission planners, and decision makers. NUOPC partners are working toward a common model architecture - a standard way of building models - in order to make it easier to collaboratively build modeling systems. To this end, they have developed the NUOPC Layer that defines conventions and a set of generic components for building coupled models using the Earth System Modeling Framework (ESMF).

SCHISM can be such a model component within an NUOPC-coupled system.  A so-called "cap" wraps SCHISM and exposes it via the NUOPC Application Programming Interface (API).  The interfaces exposed through the API are
 * import of fields
 * export of fields
 * control structure

### NUOPC - CoastalApp

The NOAA Environmental Modeling System (NEMS) Coastal Application "CoastalApp" is a NUOPC-based coupled system using the NEMS coupler.

The NUOPC cap is currently being integrated as an OCN component of the CoastalApp; this is available from a public repository [https://github.com/noaa-ocs-modeling/CoastalApp](https://github.com/noaa-ocs-modeling/CoastalApp)

#### Obtaining and building CoastalApp

```
export COASTALAPP_DIR=/my/path/to/coastalapp
git clone https://github.com/noaa-ocs-modeling/CoastalApp $COASTALAPP_DIR
cd $COASTALAPP_DIR
git checkout feature/schism
git submodule update --init --recursive SCHISM NEMS
bash ./build.sh -component "SCHISM"
```

You can add components like `WW3` or `ADCIRC`, and you may be required to choose a compiler or platform.  Consult `./build.sh -h` for help and further information.

### ESMF
The Earth System Modeling Framework (ESMF) is a suite of software tools for developing high-performance, multi-component Earth science modeling applications. Such applications may include a few or dozens of components representing atmospheric, oceanic, terrestrial, or other physical domains, and their constituent processes (dynamical, chemical, biological, etc.). Often these components are developed by different groups independently, and must be “coupled” together using software that transfers and transforms data among the components in order to form functional simulations.

SCHISM can be such a component within an ESMF-coupled system.  A so-called "cap" wraps SCHISM and exposes it via the ESMF Application Programming Interface (API).  The interfaces exposed through the API are
 * import of fields
 * export of fields
 * control structure

### ESMF - MOSSCO

MOSSCO, the "Modular System for Shelves and Coasts" is a framework for coupling
processes or domains that are originally developed in standalone numerical models.
The software MOSSCO implements this infrastructure in the form of a library of
components and couplers, and of example coupled applications. The components
"wrap" external models used in coastal and shelf sciences; these wrapped components are then coupled
to each other in the Earth System Modeling Framework (ESMF).

The [SCHISM ESMF](esmf.html) cap integrates with MOSSCO.

#### Obtaining and building MOSSCO

```
export SCHISM_BUILD_DIR=/my/path/to/schism/build
export SCHISM_ESMF_DIR=/my/path/to/schism-esmf
export MOSSCO_DIR=/my/path/to/mossco

git clone https://git.code.sf.net/p/mossco/code $MOSSCO_DIR

cd $MOSSCO_DIR
make all install

bash ./build.sh -component "SCHISM"
```

#### Using SCHISM as part of a MOSSCO coupled system

A simple preconfigured application is available in `$MOSSCO_DIR/examples/esmf/schism`. To build it, run `make` in that directory.  You can use the resulting executable as a drop-in replacement for SCHISM's standalone `pschism` executable, but you need to add the `mossco.cfg` resource file which overrides `param.nml` for control parameters of the coupled system, like start and stop time.

### ESMF - Pre-configured executables

While the API exposed through the SCHISM-ESMF library is the core functionality of the cap, there are also several
 pre-configured simple coupled systems available for you to  `make` or use as a template:

```
concurrent_esmf_test.F90
multi_schism.F90
schism_driver_interfaces.F90
schism_pdaf.F90
triple_schism.F90
```
#### Applications

(1) The SCHISM ESMF cap is used to couple SCHISM to the Parallel Data Assimilation Framework leveraging the control structures of ESMF, see `schism_pdaf.F90`.

(2) The SCHISM ESMF cap is used in the Modular System for Shelves and Coasts, see  [MOSSCO](mossco.html).  Within that system, SCHISM can be flexibly coupled to components for the atmosphere, waves, ocean BGC, generic input/output, sediment.

## Publications
Hao-Cheng Yu, Yinglong J. Zhang, Nerger Lars, Carsten Lemmen, Jason C.S Yu, Tzu-Yin Chou, Chi-Hao Chu, and Chuen-Teyr Terng: Development of a flexible data assimilation method in a 3D unstructured-grid ocean model under Earth System Modeling Framework, submitted to Geoscientific Model Development, March 2022


