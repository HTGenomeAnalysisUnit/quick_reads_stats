#Compile a static build with no dependencies on hts lib
#See: https://github.com/brentp/hts-nim#static-builds
#You need to bind your local nim packages folder as nim_pkg_folder
pkg_name="qrs"
nim_pkg_folder=$1

singularity exec \
--bind $PWD \
--bind $PWD:/load/ \
--bind $nim_pkg_folder \
docker://brentp/musl-hts-nim:latest \
/usr/local/bin/nsb \
-s $PWD/src/${pkg_name}.nim \
--nimble-file $PWD/${pkg_name}.nimble -- -d:danger -d:release
