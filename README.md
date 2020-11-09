# multi
Multi-incremental multi-resolution variational toy system

## Documentation
This code implements the final example given in the note: [Multi-incremental multi-resolution variational method: practical constraints](doc/multi.pdf)

## Getting started

- Clone ecbuild:

    cd ${path_to_source}
    git clone https://github.com/ecmwf/ecbuild

- Add ecbuild to your PATH:

   export PATH=$PATH:${path_to_source}/ecbuild/bin

- Clone multi:

   cd ${path_to_source}
   git clone https://github.com/benjaminmenetrier/multi

- Compile the code:

   cd ${path_to_build}
   ecbuild ${path_to_source}/multi
   make

- Run the code:

   ${path_to_build}/bin/multi < ${path_to_source}/namelist
