PYTHON = python3

SCRIPTS_DIR = hadiahmad
SCRIPTS = $(SCRIPTS_DIR)/DataProcessing.py $(SCRIPTS_DIR)/MyHPTest.py $(SCRIPTS_DIR)/MyBonferri.py $(SCRIPTS_DIR)/MyFDR.py
OUTPUT_DIR = output
OUTPUTS = $(OUTPUT_DIR)/ArrayData.txt $(OUTPUT_DIR)/SeqData.txt

1:
	$(PYTHON) hadiahmad/DataProcessing.py
	@echo "---------------\nDataProcessing.py has been run\n---------------"

2:
	for output in $(OUTPUTS); do $(PYTHON) $(SCRIPTS_DIR)/MyHPTest.py $$output; done
	@echo "---------------\nMyHPTest.py has been run\n---------------"

#runs both
# python3 hadiahmad/MyBonferri.py output/ArrayData.txt
# python3 hadiahmad/MyBonferri.py output/SeqData.txt
3:
	for output in $(OUTPUTS); do $(PYTHON) $(SCRIPTS_DIR)/MyBonferri.py $$output; done
	@echo "---------------\nMyBonferri.py has been run\n---------------"

# runs both
# python3 hadiahmad/MyFDR.py output/ArrayData.txt 20
# python3 hadiahmad/MyFDR.py output/SeqData.txt 20
4:
	for output in $(OUTPUTS); do $(PYTHON) $(SCRIPTS_DIR)/MyFDR.py $$output 20; done
	for output in $(OUTPUTS); do $(PYTHON) $(SCRIPTS_DIR)/MyFDR.py $$output 50; done
	for output in $(OUTPUTS); do $(PYTHON) $(SCRIPTS_DIR)/MyFDR.py $$output 100; done
	for output in $(OUTPUTS); do $(PYTHON) $(SCRIPTS_DIR)/MyFDR.py $$output 200; done
	@echo "---------------\nMyFDR.py has been run\n---------------"

all: 1 2 3 4
	@echo "All .py files have been run"

%.py:
	$(PYTHON) $@

clean:
	rm -f *.pyc