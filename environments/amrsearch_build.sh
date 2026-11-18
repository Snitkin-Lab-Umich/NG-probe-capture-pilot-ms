#!/usr/bin/env bash
set -eo pipefail

cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/environments

module load mamba

source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^amrsearch "; then
  mamba env create -f amrsearch.yaml
fi
conda activate amrsearch

git --version
mvn -version
java -version
makeblastdb -version

cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/containers

# re resolve for if the library ever needs updating
AMRSEARCH_COMMIT=767320f0c0d63a8184e7c9d58aee00b37f31d380
AMRLIBRARIES_COMMIT=f8806f3c3d4fa49fb5bfbfd86898a38a3515cfd4

mkdir -p amr-search
cd amr-search
git init -q
git remote add origin https://github.com/pathogenwatch-oss/amr-search
git fetch --depth=1 origin "$AMRSEARCH_COMMIT"
git checkout -q FETCH_HEAD

mkdir -p libraries/amr-libraries
cd libraries/amr-libraries
git init -q
git remote add origin https://github.com/pathogenwatch/amr-libraries
git fetch --depth=1 origin "$AMRLIBRARIES_COMMIT"
git checkout -q FETCH_HEAD
cd ../..

mkdir -p libraries/amr-test-libraries
cp libraries/amr-libraries/antimicrobials.toml libraries/amr-test-libraries/antimicrobials.toml

mvn -s /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/environments/maven-settings.xml -Dlicense.skipDownloadLicenses=true install

cp resources/taxid.map build/databases/taxid.map

cat > BUILD_INFO.txt << EOF
amr-search commit: $AMRSEARCH_COMMIT
amr-libraries commit: $AMRLIBRARIES_COMMIT
build date: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

# Notes:
# - Run this whole script via srun on a compute node; the login node kills conda/mvn as heavy processes
# - pom.xml actually needs Java 11, despite the README saying Java 8
# - maven-settings.xml supplies the Great Lakes proxy; mvn ignores http_proxy/https_proxy env vars
# - amr-test-libraries is a private submodule we can't access; the placeholder dir+file satisfies the builder's existence check without contributing any data
# - license.skipDownloadLicenses avoids multi-minute hangs on dead license URLs during the build
# - taxid.map is generated under resources/ but the runner expects it in build/databases/, hence the copy step
# - commits are pinned (not -b main) so rebuilds are reproducible; BUILD_INFO.txt records what was actually built
