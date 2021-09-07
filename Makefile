TWL=test-small.wdl
BUNDLE=bundle.zip
CROMWELL_JAR=~/bin/cromwell-54.jar
WOMTOOL_JAR=~/bin/womtool-54.jar

METAG_FULL_WORKFLOW=nmdc-metag_full.wdl
METAG_INPUT=inputs.json
METAG_LABELS=labels.json

METAT_FULL_WORKFLOW=nmdc-metat_full.wdl
METAT_INPUT=inputs_metat.json
METAT_LABELS=labels_metat.json

zip:
	zip $(BUNDLE) *wdl

validate:
	java -jar $(WOMTOOL_JAR) validate $(TWL)

test:
	echo "use submit"

submit:
	java -jar $(CROMWELL_JAR) submit -h $(CROMWELL_URL) $(TWL) -p $(BUNDLE)

validatemeta:
	java -jar $(WOMTOOL_JAR) validate $(METAG_FULL_WORKFLOW) -i $(METAG_INPUT)

metag:
	java -jar $(CROMWELL_JAR) submit -h $(CROMWELL_URL) $(METAG_FULL_WORKFLOW) -p $(BUNDLE) -i $(METAG_INPUT) -l $(METAG_LABELS)

validatemetat:
	java -jar $(WOMTOOL_JAR) validate $(METAT_FULL_WORKFLOW) -i $(METAT_INPUT)

metat:
	java -jar $(CROMWELL_JAR) submit -h $(CROMWELL_URL) $(METAT_FULL_WORKFLOW) -p $(BUNDLE) -i $(METAT_INPUT) -l $(METAT_LABELS)

