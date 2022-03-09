.PHONY: html

HTMLPATH?=$(pwd)./notebooks/

# initialize virtual environment in local folder
# same as:
# $ virtualenv -p /usr/local/bin/python3 venv3.37
# $ source venv3/bin/activate
# $ python3 -m pip install -r requirements.txt
init:
	( \
		virtualenv -p /usr/bin/python3 venv3; \
  	. ./venv3/bin/activate; \
		python --version; \
		python3 -m pip install -r requirements.txt; \
	)

# install missing packages newly added to requirements.txt
install:
	( \
    . ./venv3/bin/activate; \
    python --version; \
		python3 -m pip install --upgrade pip; \
    python3 -m pip install -r requirements.txt; \
  )

# start jupyter notebook in virtual environment 
start:
	( \
		. venv3/bin/activate; \
		jupyter nbextension enable --py --sys-prefix qgrid; \
		jupyter nbextension enable --py --sys-prefix widgetsnbextension; \
		jupyter notebook \
	)

# Use as follows:
# $ make html HTMLPATH=./value_vs_rewards/
html:
	( \
    . ./venv3/bin/activate; \
		jupyter nbconvert --to html $(HTMLPATH)*.ipynb \
	) 

# Extract source code form notebook to python file
# $ make extract
extract:
	( \
	 . ./venv3/bin/activate; \
	NBFILE="./notebooks/block_profit.ipynb" PYFILE="./src/profitability/profitability.py" python3 ./src/util/nb_extract.py \
	)	 

# Generate the tables 
tables: tblexlclassic tblexlmarkov

tblexlclassic:
	( \
	 . ./venv3/bin/activate; \
	python3 ./src/visualization/tbl_profits_exclusion_classic.py > ./paper/preprint/texsrc/tbl_profits_exclusion_classic.tex; \
	)

tblexlmarkov:
	( \
	 . ./venv3/bin/activate; \
	python3 ./src/visualization/tbl_profits_exclusion_markov.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_0.tex; \
	P_V="0.1" KR="3" python3 ./src/visualization/tbl_profits_exclusion_markov.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_pv10.tex; \
	P_V="0.1" KR="3" N_FORK="1" python3 ./src/visualization/tbl_profits_exclusion_markov.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_pv10nf1.tex; \
	P_V="0.1" KR="3" N_MAIN="1" python3 ./src/visualization/tbl_profits_exclusion_markov.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_pv10nm1.tex; \
	)

tblexlmarkoveffort:
	( \
	 . ./venv3/bin/activate; \
	P_V="0.1" KR="3" N_FORK="1" python3 ./src/visualization/tbl_profits_exclusion_markov_effort.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_effort_pv10nf1.tex; \
	P_V="0.1" KR="3" N_MAIN="1" python3 ./src/visualization/tbl_profits_exclusion_markov_effort.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_effort_pv10nm1.tex; \
	P_V="0.1" KR="3" python3 ./src/visualization/tbl_profits_exclusion_markov_effort.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_effort_pv10.tex; \
	P_V="0.1" KR="3" Z="6" python3 ./src/visualization/tbl_profits_exclusion_markov_effort.py > ./paper/preprint/texsrc/tbl_profits_exclusion_markov_effort_pv10z6.tex; \
	)

clean:
	-rm -i $(HTMLPATH)*.html

