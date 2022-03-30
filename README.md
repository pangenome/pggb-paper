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
latexmk -pdf sn-article.tex
```
