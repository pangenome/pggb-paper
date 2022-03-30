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

## submission guidelines

From [here](https://www.nature.com/nbt/content):

A Brief Communication reports a concise study of high quality and broad interest.

Format
- Brief unreferenced abstract – 3 sentences, up to 70 words, which will appear on Medline.
- Title – up to 10 words (or 90 characters).
- Main text – 1,000-1,500 words, including abstract, references and figure legends, and contains no headings.
- Display items – up to 2 items, although this may be flexible at the discretion of the editor, provided the page limit is observed.
- Online Methods section should be included
- References – as a guideline, we typically recommend up to 20. Article titles are omitted from the reference list.
- Brief Communications should include received/accepted dates.
- Brief Communications may be accompanied by supplementary information.
- Brief Communications are peer reviewed.
