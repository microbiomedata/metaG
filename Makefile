BUNDLE=bundle.zip
CROMWELL_JAR=~/bin/cromwell-54.jar
WOMTOOL_JAR=~/bin/womtool-54.jar

METAG_FULL_WORKFLOW=nmdc-metag_full.wdl
METAG_INPUT=inputs.json
METAG_LABELS=labels.json

zip:
	(cd subworkflows && zip ../$(BUNDLE) *wdl)

prep:
	mkdir -p subworkflows
	git clone --depth 1 https://github.com/microbiomedata/ReadsQC -b b1.0.3 subworkflows/ReadsQC
	git clone --depth 1 https://github.com/microbiomedata/metaAssembly -b b1.0.3 subworkflows/metaAssembly
	git clone --depth 1 https://github.com/microbiomedata/mg_annotation -b add_outputs subworkflows/mg_annotation
	git clone --depth 1 https://github.com/microbiomedata/ReadbasedAnalysis -b b1.0.3 subworkflows/ReadbasedAnalysis
	git clone --depth 1 https://github.com/microbiomedata/metaMAGs -b b1.0.3 subworkflows/metaMAGs
	(cd subworkflows && ln -s */*wdl .)

validate:
	java -jar $(WOMTOOL_JAR) validate $(METAG_FULL_WORKFLOW) -i $(METAG_INPUT)

test:
	java -jar $(CROMWELL_JAR) submit -h $(CROMWELL_URL) $(METAG_FULL_WORKFLOW) -p $(BUNDLE) -i $(METAG_INPUT) -l $(METAG_LABELS)

testarea:
	rm -rf ./t
	mkdir ./t
	(cd t && ln -s ../subworkflows/*.wdl ../*.wdl ../tests/*.wdl .)

testmetadata:
	java -jar /tmp/cromwell.jar run ./t/validate_finish.wdl
