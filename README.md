# pggb-paper

## building the manuscript

```shell
# Dependencies
sudo apt-get -y install texlive texlive-latex-recommended \
        texlive-pictures texlive-latex-extra texlive-fonts-extra texlive-font-utils\
        texlive-science

sudo apt install latexmk

git clone https://github.com/pangenome/pggb-paper
cd pggb-paper/manuscript
#latexmk -pdf sn-article.tex
pdflatex article && bibtex article && pdflatex article && pdflatex article
```

## submission guidelines

From [here](https://www.nature.com/nmeth/content):

### Brief Communication
A Brief Communication is a concise report describing potentially groundbreaking yet preliminary method or tool developments, highly practical tweaks to an existing method or tool, software platforms, resources of broad interest, and technical critiques of widely used methodologies.

#### Format

Abstract – up to 70 words, unreferenced.
Main text – 1,200 words (up to 1600 words with editorial discretion), including abstract, references and figure legends.
The main text should not contain sections nor subheadings. 
Display items – maximum 2 figures and/or tables (up to 3 items with editorial discretion).
Online Methods section should be included and should contain subheadings.
References – as a guideline, we typically recommend up to 20. 
Brief Communications include received/accepted dates. 
Brief Communications may be accompanied by supplementary information. 
Brief Communications are peer reviewed.
