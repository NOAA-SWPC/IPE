bucket = https://noaa-nws-wam-ipe-pds.s3.amazonaws.com

data:
	wget $bucket/NON-OPERATIONAL_RESEARCH/SWORD_DATA/IPE_STANDALONE_FIX.tar.gz
	tar xvf IPE_STANDALONE_FIX.tar.gz
	rm -rf  IPE_STANDALONE_FIX.tar.gz

regression:
	wget $bucket/NON-OPERATIONAL_RESEARCH/SWORD_DATA/IPE_STANDALONE_REGRESSION.tar.gz
	tar xvf IPE_STANDALONE_REGRESSION.tar.gz
	rm -rf  IPE_STANDALONE_REGRESSION.tar.gz

NP=8

help:
	@echo "make -j test          - run full test on $(NP) processors"
	@echo "make -j test NP=4     - run test on 4 processors"
	@echo "make -j test_compile  - (re)compile ipe.x"
	@echo "make test_rundir      - create run directory 'run/'"
	@echo "make test_run NP=3    - cd run; mpiexec -n 3 ipe.x"
	@echo "make test_check       - check test results"
	@echo "make clean            - remove files created by test"

test:
	@echo test_compile > test.diff
	$(MAKE) test_compile
	@echo test_rundir >> test.diff
	$(MAKE) test_rundir
	@echo test_run >> test.diff
	$(MAKE) test_run
	@echo test_check >> test.diff
	$(MAKE) test_check

test_compile:
	cd ..; make -j; make install

test_rundir: data
	rm -rf run
	mkdir -p run
	cd run; ln -s ../data/* .
	cd run; ln -s ../../install/bin/ipe.x .

TESTFILE = IPE_State.apex.201303160005.nc

test_run:
	rm -f run/$(TESTFILE)
	cd run;	mpiexec -n $(NP) ./ipe.x > runlog

test_check: regression
	scripts/compare.py run/$(TESTFILE) regression/$(TESTFILE).derecho \
	> test.diff
	@cat test.diff

clean:
	rm -rf run test.diff
