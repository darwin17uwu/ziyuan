# Compila la presentación Beamer del Capítulo 4:
# "Regresión por Mínimos Cuadrados" (Brunton, Optimization Bootcamp).
#
# Uso:
#   make          # compila el PDF
#   make clean    # borra los archivos auxiliares
#   make distclean# borra auxiliares y PDF
#   make rebuild  # clean + distclean + compila de nuevo
#   make check    # muestra el PDF generado
#   make preview  # renderiza cada página a PNG en preview/
#
# Requiere: pdflatex (TeX Live / MiKTeX / MacTeX) y Poppler para preview.

TEX := ch4_minimos_cuadrados.tex
PDF := $(TEX:.tex=.pdf)
LOG := $(TEX:.tex=.log)

PASSES := 5

.PHONY: all clean distclean rebuild check preview

all: $(PDF)

# pdflatex se repite hasta que el índice y los marcadores se estabilicen
$(PDF): $(TEX)
	@echo "Compilando $(PDF)..."
	@i=1; \
	while [ $$i -le $(PASSES) ]; do \
	  echo "  pdflatex, pasada $$i/$(PASSES)"; \
	  pdflatex -interaction=nonstopmode -halt-on-error $(TEX) >/dev/null || exit 1; \
	  if ! grep -q "Rerun to get" $(LOG); then \
	    echo "  índice estable tras $$i pasada(s)"; \
	    break; \
	  fi; \
	  i=$$((i + 1)); \
	done
	@echo "Listo: $(PDF)"

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.fls *.fdb_latexmk *.synctex.gz

distclean: clean
	rm -f $(PDF)

rebuild: distclean all

check: $(PDF)
	@ls -lh $(PDF)
	@pdfinfo $(PDF) 2>/dev/null | grep -i '^Pages' || true

preview: $(PDF)
	mkdir -p preview
	@for page in $$(seq 1 $$(pdfinfo $(PDF) | awk '/^Pages:/ {print $$2}')); do \
	  pdftoppm -f $$page -l $$page -scale-to 1400 -singlefile -png $(PDF) preview/slide-$$page; \
	done
	@echo "Preview: $$(ls preview/*.png | wc -l) paginas en preview/"
