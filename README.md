# IDPs.hydro
####This project could not have been possible without the collaboration of [Swarnadeep Seth](targetURL "https://www.linkedin.com/in/swarnadeepseth/") and [their extensive knowledge in related fields](targetURL "https://swarnadeepseth.github.io/")
This project has two main sections; Phase 1 and Phase 2.
##Phase 1
This is the data gathering phase, all the data is calculated and organized for the user.
Note: Phase 1 is designed to be run on a super computer using slurm
Phase 1 does require two unique environments to be created to ensure working capabilies; The environments are hoomd_w_panda and my_env.
If you do not care for these names or need to specify your particular modules, you will need to change the run_code.py file commands list(located 159-165). 
hoomd_w_panda is for organization of data and data computation whilst my_env is for analysis of data to produce information more understandable.
Lines of interest for doing so for hoomd_w_panda are 160,161, and 164.
###hoomd_w_panda packages,version,build and Channel(if no Channel specified then should be downloadable with pip)
 Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main
attrs                     19.3.0                     py_0
backcall                  0.1.0                    py37_0
blas                      1.0                         mkl
bleach                    3.1.0                    py37_0
ca-certificates           2022.07.19           h06a4308_0
certifi                   2022.9.14        py37h06a4308_0
cudatoolkit               10.1.168                      0
cycler                    0.10.0                   py37_0
dbus                      1.13.12              h746ee38_0
decorator                 4.4.1                      py_0
defusedxml                0.6.0                      py_0
embree3                   3.6.1                hc9558a2_0    conda-forge
entrypoints               0.3                      py37_0
expat                     2.2.6                he6710b0_0
fontconfig                2.13.0               h9420a91_0
freetype                  2.9.1                h8a8886c_1
fresnel                   0.11.0           py37h9de70de_0    conda-forge
freud                     2.0.0            py37hc9558a2_0    conda-forge
glib                      2.56.2               hd408876_0
gmp                       6.1.2                h6c8ec71_1
gsd                       1.9.3            py37hc1659b7_0    conda-forge
gst-plugins-base          1.14.0               hbbd80ab_1
gstreamer                 1.14.0               hb453b48_1
hoomd                     2.8.0            py37h9de70de_0    conda-forge
icu                       58.2                 h9c2bf20_1
importlib_metadata        0.23                     py37_0
intel-openmp              2019.4                      243
ipykernel                 5.1.3            py37h39e3cac_0
ipython                   7.9.0            py37h39e3cac_0
ipython_genutils          0.2.0                    py37_0
ipywidgets                7.5.1                      py_0
jedi                      0.15.1                   py37_0
jinja2                    2.10.3                     py_0
jpeg                      9b                   h024ee3a_2
jsonschema                3.1.1                    py37_0
jupyter                   1.0.0                    py37_7
jupyter_client            5.3.4                    py37_0
jupyter_console           6.0.0                    py37_0
jupyter_core              4.6.0                    py37_0
kiwisolver                1.1.0            py37he6710b0_0
libedit                   3.1.20181209         hc058e9b_0
libffi                    3.2.1                hd88cf55_4
libgcc-ng                 9.1.0                hdf63c60_0
libgfortran-ng            7.3.0                hdf63c60_0
libpng                    1.6.37               hbc83047_0
libsodium                 1.0.16               h1bed415_0
libstdcxx-ng              9.1.0                hdf63c60_0
libtiff                   4.0.10               h2733197_2
libuuid                   1.0.3                h1bed415_2
libxcb                    1.13                 h1bed415_1
libxml2                   2.9.9                hea5a465_1
markupsafe                1.1.1            py37h7b6447c_0
matplotlib                3.1.1            py37h5429711_0
mistune                   0.8.4            py37h7b6447c_0
mkl                       2019.4                      243
mkl-service               2.3.0            py37he904b0f_0
mkl_fft                   1.0.14           py37ha843d7b_0
mkl_random                1.1.0            py37hd6b4f25_0
more-itertools            7.2.0                    py37_0
mpi                       1.0                     openmpi    conda-forge
nbconvert                 5.6.0                    py37_1
nbformat                  4.4.0                    py37_0
ncurses                   6.1                  he6710b0_1
notebook                  6.0.1                    py37_0
numpy                     1.17.2           py37haad9e8e_0
numpy-base                1.17.2           py37hde5b4d6_0
olefile                   0.46                     py37_0
openmpi                   4.0.1                hc99cbb1_2    conda-forge
openssl                   1.1.1q               h7f8727e_0
pandas                    0.25.1                   pypi_0    pypi
pandoc                    2.2.3.2                       0
pandocfilters             1.4.2                    py37_1
parso                     0.5.1                      py_0
pcre                      8.43                 he6710b0_0
pexpect                   4.7.0                    py37_0
pickleshare               0.7.5                    py37_0
pillow                    6.2.0            py37h34e0f95_0
pip                       19.3.1                   py37_0
prometheus_client         0.7.1                      py_0
prompt_toolkit            2.0.10                     py_0
ptyprocess                0.6.0                    py37_0
pygments                  2.4.2                      py_0
pyparsing                 2.4.2                      py_0
pyqt                      5.9.2            py37h05f1152_2
pyrsistent                0.15.4           py37h7b6447c_0
python                    3.7.5                h0371630_0
python-dateutil           2.8.0                    py37_0
pytz                      2019.3                     py_0
pyzmq                     18.1.0           py37he6710b0_0
qhull                     2015.2            h6bb024c_1000    conda-forge
qt                        5.9.7                h5867ecd_1
qtconsole                 4.5.5                      py_0
readline                  7.0                  h7b6447c_5
rowan                     1.2.2                      py_1    conda-forge
scipy                     1.3.1            py37h7c811a0_0
send2trash                1.5.0                    py37_0
setuptools                41.6.0                   py37_0
sip                       4.19.8           py37hf484d3e_0
six                       1.12.0                   py37_0
sqlite                    3.30.1               h7b6447c_0
tbb                       2019.9               hc9558a2_0    conda-forge
tbb-devel                 2019.9               hc9558a2_0    conda-forge
terminado                 0.8.2                    py37_0
testpath                  0.4.2                    py37_0
tk                        8.6.8                hbc83047_0
tornado                   6.0.3            py37h7b6447c_0
traitlets                 4.3.3                    py37_0
wcwidth                   0.1.7                    py37_0
webencodings              0.5.1                    py37_1
wheel                     0.33.6                   py37_0
widgetsnbextension        3.5.1                    py37_0
xz                        5.2.4                h14c3975_4
zeromq                    4.3.1                he6710b0_3
zipp                      0.6.0                      py_0
zlib                      1.2.11               h7b6447c_3
zstd                      1.3.7                h0b5b093_0




Lines of interest for my_env are 162 and 163.
###my_envc packages,version(does not have conda list output)

absl-py                            0.9.0
alabaster                          0.7.12
anaconda-client                    1.7.2
anaconda-navigator                 1.9.12
anaconda-project                   0.8.3
argh                               0.26.2
argon2-cffi                        20.1.0
asn1crypto                         1.4.0
astroid                            2.4.2
astropy                            4.0.1.post1
async-generator                    1.10
atomicwrites                       1.4.0
attrs                              20.2.0
autopep8                           1.5.4
Babel                              2.8.0
backcall                           0.2.0
backports.functools-lru-cache      1.6.1
backports.shutil-get-terminal-size 1.0.0
backports.tempfile                 1.0
backports.weakref                  1.0.post1
basemap                            1.2.1
beautifulsoup4                     4.9.3
biopython                          1.81
bitarray                           1.5.3
bkcharts                           0.2
bleach                             3.2.1
bokeh                              2.2.1
boto                               2.49.0
Bottleneck                         1.3.2
brotlipy                           0.7.0
certifi                            2022.5.18.1
cffi                               1.14.3
chardet                            3.0.4
click                              7.1.2
cloudpickle                        1.6.0
clyent                             1.2.2
colorama                           0.4.3
conda                              4.13.0
conda-build                        3.18.11
conda-package-handling             1.7.0
conda-verify                       3.4.2
contextlib2                        0.6.0.post1
cryptography                       3.1.1
cycler                             0.10.0
Cython                             0.29.21
cytoolz                            0.11.0
dask                               2.30.0
DateTime                           5.1
decorator                          4.4.2
defusedxml                         0.6.0
diff-match-patch                   20200713
distributed                        2.30.0
docutils                           0.16
entrypoints                        0.3
et-xmlfile                         1.0.1
fastcache                          1.1.0
fasteners                          0.18
filelock                           3.0.12
flake8                             3.8.4
Flask                              1.1.2
fsspec                             0.8.0
future                             0.18.2
gevent                             20.9.0
glob2                              0.7
gmpy2                              2.0.8
greenlet                           0.4.17
GridDataFormats                    1.0.1
gsd                                2.9.0
h5py                               2.10.0
HeapDict                           1.0.1
html5lib                           1.1
idna                               2.10
imageio                            2.9.0
imagesize                          1.2.0
importlib-metadata                 2.0.0
iniconfig                          0.0.0
intervaltree                       3.1.0
ipykernel                          5.3.4
ipython                            7.18.1
ipython-genutils                   0.2.0
ipywidgets                         7.5.1
isort                              5.6.4
itsdangerous                       1.1.0
jdcal                              1.4.1
jedi                               0.17.1
jeepney                            0.4.3
Jinja2                             2.11.2
joblib                             0.17.0
json5                              0.9.5
jsonschema                         3.2.0
jupyter                            1.0.0
jupyter-client                     6.1.7
jupyter-console                    6.2.0
jupyter-core                       4.6.3
jupyterlab                         2.2.6
jupyterlab-pygments                0.1.2
jupyterlab-server                  1.2.0
keyring                            21.4.0
kiwisolver                         1.2.0
lazy-object-proxy                  1.4.3
libarchive-c                       2.9
llvmlite                           0.34.0
locket                             0.2.0
lxml                               4.5.2
mamba                              0.15.3
MarkupSafe                         1.1.1
matplotlib                         3.3.1
mccabe                             0.6.1
MDAnalysis                         2.2.0
mimesis                            4.0.0
mistune                            0.8.4
mkl-fft                            1.2.0
mkl-random                         1.1.1
mkl-service                        2.3.0
mmtf-python                        1.1.3
mock                               4.0.2
more-itertools                     8.5.0
mpmath                             1.1.0
mrcfile                            1.4.3
msgpack                            1.0.0
multipledispatch                   0.6.0
navigator-updater                  0.2.1
nbclient                           0.5.0
nbconvert                          6.0.7
nbformat                           5.0.7
nest-asyncio                       1.4.1
networkx                           2.5
nltk                               3.5
nose                               1.3.7
notebook                           6.1.4
numba                              0.51.2
numexpr                            2.7.1
numpy                              1.19.1
numpydoc                           1.1.0
olefile                            0.46
openpyxl                           3.0.5
packaging                          20.4
pandas                             1.1.3
pandocfilters                      1.4.2
parso                              0.7.0
partd                              1.1.0
path                               15.0.0
pathlib2                           2.3.5
pathtools                          0.1.2
patsy                              0.5.1
pep8                               1.7.1
pexpect                            4.8.0
pickleshare                        0.7.5
Pillow                             7.2.0
pip                                20.2.3
pkginfo                            1.5.0.1
pluggy                             0.13.1
ply                                3.11
prometheus-client                  0.8.0
prompt-toolkit                     3.0.7
psutil                             5.7.2
ptyprocess                         0.6.0
py                                 1.9.0
pycodestyle                        2.6.0
pycosat                            0.6.3
pycparser                          2.20
pycurl                             7.43.0.6
pydocstyle                         5.1.1
pyflakes                           2.2.0
Pygments                           2.7.1
pylint                             2.6.0
pyodbc                             4.0.0-unsupported
pyOpenSSL                          19.1.0
pyparsing                          2.4.7
pyproj                             2.6.1.post1
pyrsistent                         0.17.3
pyshp                              2.1.2
PySocks                            1.7.1
pytest                             0.0.0
python-dateutil                    2.8.1
python-jsonrpc-server              0.4.0
python-language-server             0.35.1
pytorch-triton                     2.1.0+440fd1bf20
pytz                               2020.1
PyWavelets                         1.1.1
pyxdg                              0.26
PyYAML                             5.3.1
pyzmq                              19.0.2
QDarkStyle                         2.8.1
QtAwesome                          1.0.1
qtconsole                          4.7.7
QtPy                               1.9.0
regex                              2020.9.27
requests                           2.24.0
rope                               0.18.0
Rtree                              0.9.4
ruamel-yaml                        0.15.87
scikit-image                       0.16.2
scikit-learn                       0.23.2
scipy                              1.5.2
seaborn                            0.11.0
SecretStorage                      3.1.2
Send2Trash                         1.5.0
setuptools                         50.3.0.post20201006
simplegeneric                      0.8.1
singledispatch                     3.4.0.3
sip                                4.19.13
six                                1.15.0
snowballstemmer                    2.0.0
sortedcollections                  1.2.1
sortedcontainers                   2.2.2
soupsieve                          2.0.1
Sphinx                             3.2.1
sphinxcontrib-applehelp            1.0.2
sphinxcontrib-devhelp              1.0.2
sphinxcontrib-htmlhelp             1.0.3
sphinxcontrib-jsmath               1.0.1
sphinxcontrib-qthelp               1.0.3
sphinxcontrib-serializinghtml      1.1.4
sphinxcontrib-websupport           1.2.4
spyder                             4.1.5
spyder-kernels                     1.9.4
SQLAlchemy                         1.3.19
statsmodels                        0.12.0
sympy                              1.6.2
tables                             3.6.1
tblib                              1.7.0
terminado                          0.9.1
testpath                           0.4.4
threadpoolctl                      2.1.0
toml                               0.10.1
toolz                              0.11.1
torch                              2.1.0.dev20230623+cu121
torchaudio                         2.1.0.dev20230623+cu121
torchfile                          0.1.0
torchvision                        0.16.0.dev20230623+cu121
tornado                            6.0.4
tqdm                               4.50.2
traitlets                          5.0.4
typing-extensions                  3.7.4.3
ujson                              4.0.1
unicodecsv                         0.14.1
urllib3                            1.25.10
visdom                             0.1.8.9
watchdog                           0.10.3
wcwidth                            0.2.5
webencodings                       0.5.1
websocket-client                   0.58.0
Werkzeug                           1.0.1
wheel                              0.35.1
widgetsnbextension                 3.5.1
wrapt                              1.11.2
wurlitzer                          2.0.1
xlrd                               1.2.0
XlsxWriter                         1.3.6
xlwt                               1.3.0
xmltodict                          0.12.0
yapf                               0.30.0
zict                               2.0.0
zipp                               3.3.0
zope.event                         4.4
zope.interface                     5.1.2
