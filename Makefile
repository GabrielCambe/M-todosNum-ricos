NOME = retroSub
SRC = ./src/
INCLUDE = ./include/
MAIN_src = $(SRC)main.c

all:
	rm $(NOME) || true
bck_sub: all
	gcc -o $(NOME) retroSubstituicao.c -Wall -O0 -D BACK_SUB
fwd_sub: all
	gcc -o $(NOME) retroSubstituicao.c -Wall -O0 -D FOWR_SUB
purge:
	rm *~ $(NOME) || true

#regras para o controle de versão
commit:
	git update-index
	git commit -m "commit com MAKE"
	git push

update:
	git update-index --add
	git commit -a --dry-run
	git add --all
