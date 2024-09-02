help([[
  loads GSMWAM-IPE prerequisites on WCOSS2
]])

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.1.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

craype_ver=os.getenv("craype_ver") or "2.7.8"
load(pathJoin("craype", craype_ver))

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.7"
load(pathJoin("cray-mpich", cray_mpich_ver))

cray_libsci_ver=os.getenv("cray_libsci_ver") or "21.06.1.1"
load(pathJoin("cray-libsci", cray_libsci_ver))

cray_pals_ver=os.getenv("cray_pals_ver") or "1.0.12"
load(pathJoin("cray-pals", cray_pals_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("netcdf", netcdf_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.2"
load(pathJoin("nemsio", nemsio_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.10.6"
load(pathJoin("hdf5", hdf5_ver))

pnetcdf_ver=os.getenv("pnetcdf_ver") or "1.12.2"
load(pathJoin("pnetcdf", pnetcdf_ver))

sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))

pio_ver=os.getenv("pio_ver") or "2.5.10"
load(pathJoin("pio", pio_ver))

esmf_ver=os.getenv("esmf_ver") or "8.4.1"
load(pathJoin("esmf", esmf_ver))

-- comio load occurs out of $HOMEwfs/modulefiles and is implicitly specified
-- by build.ver as the version of COMIO that is checked out comes from build.ver.
if os.getenv("HOMEwfs") then
  prepend_path("MODULEPATH", pathJoin(os.getenv("HOMEwfs"), "modulefiles"))
else
  append_path("MODULEPATH", "/lfs/h1/swpc/wam/noscrub/swpc.wam/wam-ipe_workflow/modulefiles")
end
comio_ver=os.getenv("comio_ver") or "v0.0.10"
load(pathJoin("comio", comio_ver))

python_ver=os.getenv("python_ver") or "3.8.6"
load(pathJoin("python", python_ver))

-- set recommended Fortran compiler flags
setenv("FCFLAGS", "-O2 -fp-model precise -ftz -fast-transcendentals -no-prec-div -no-prec-sqrt -align array64byte -align sequence")

whatis("Description: WAM-IPE build environment")

