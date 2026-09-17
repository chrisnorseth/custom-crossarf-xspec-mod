#!/usr/bin/env bash
set -euo pipefail

ENV_PREFIX="/Users/christiannorseth/yt-conda/henv"
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_NAME="${XSPEC_XCRSARF_PACKAGE:-xcpkg}"
BUILD_DIR="${SRC_DIR}/build/${PKG_NAME}"
rm -rf "${BUILD_DIR}"
TOOLBIN="${BUILD_DIR}/toolbin"

mkdir -p "${TOOLBIN}" "${BUILD_DIR}/.home"
cp "${SRC_DIR}/XCrossArf.cxx" "${BUILD_DIR}/"
cp "${SRC_DIR}/XCrossArf.h" "${BUILD_DIR}/"
cp "${SRC_DIR}/xcrsarfWrapper.cxx" "${BUILD_DIR}/"
cp "${SRC_DIR}/xcrsarfWrapper.h" "${BUILD_DIR}/"
cp "${SRC_DIR}/XCrossArfSetup.cxx" "${BUILD_DIR}/"
cp "${SRC_DIR}/lmodel_xcrsarf.dat" "${BUILD_DIR}/lmodel.dat"

cat > "${TOOLBIN}/arm64-apple-darwin20.0.0-clang" <<'WRAP'
#!/usr/bin/env bash
exec /usr/bin/clang "$@"
WRAP

cat > "${TOOLBIN}/arm64-apple-darwin20.0.0-clang++" <<'WRAP'
#!/usr/bin/env bash
exec /usr/bin/clang++ "$@"
WRAP

cat > "${TOOLBIN}/clang++" <<'WRAP'
#!/usr/bin/env bash
exec /usr/bin/clang++ "$@"
WRAP

chmod +x "${TOOLBIN}/arm64-apple-darwin20.0.0-clang" \
  "${TOOLBIN}/arm64-apple-darwin20.0.0-clang++" \
  "${TOOLBIN}/clang++"

export CONDA_PREFIX="${ENV_PREFIX}"
export ENV_PREFIX
export HOME="${BUILD_DIR}/.home"
export PATH="${TOOLBIN}:${ENV_PREFIX}/bin:${PATH}"

# shellcheck source=/dev/null
source "${ENV_PREFIX}/bin/heainit.sh"

cd "${BUILD_DIR}"
initpackage "${PKG_NAME}" lmodel.dat .

python - <<'PY'
import os
from pathlib import Path
p = Path("Makefile")
text = p.read_text()
env_prefix = os.environ["ENV_PREFIX"]
text = text.replace("${F77LIBS4C}", f"-L{env_prefix}/lib -lgfortran -lquadmath -lm")
text = text.replace("-lXSModel", "-lXSModel -lhdsp_6.36")
p.write_text(text)

init = Path("lpack_xcpkg.cxx")
text = init.read_text()
text = text.replace(
    'extern "C" int Xcpkg_SafeInit(Tcl_Interp* tclInterp);',
    'extern "C" int Xcpkg_SafeInit(Tcl_Interp* tclInterp);\n'
    'extern "C" int XCrossArfSetupCmd(ClientData, Tcl_Interp*, int, Tcl_Obj* CONST []);\n'
    'extern "C" int XCrossArfRunCmd(ClientData, Tcl_Interp*, int, Tcl_Obj* CONST []);')
text = text.replace(
    '        createxcpkgFunctionMap();',
    '        createxcpkgFunctionMap();\n'
    '        Tcl_CreateObjCommand(tclInterp, "xcrsarfsetup", '
    'XCrossArfSetupCmd, 0, 0);\n'
    '        Tcl_CreateObjCommand(tclInterp, "xcrsarfrun", '
    'XCrossArfRunCmd, 0, 0);')
init.write_text(text)
PY

hmake

cat <<EOF

Built XSPEC crossarf clone package:
  ${BUILD_DIR}/lib${PKG_NAME}.dylib

Load it in XSPEC with:
  lmod ${PKG_NAME} ${BUILD_DIR}
  model xcrsarf*<your source model>
EOF
