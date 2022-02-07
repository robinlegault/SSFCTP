current_dir = $(shell pwd)
main: main.o combo.o ssfctp.o util.o generator.o
	gcc main.o combo.o ssfctp.o util.o generator.o -o main
combo.o: combo.c combo.h
	gcc -c combo.c
generator.o: generator.c generator.h
	gcc -c generator.c
util.o: util.c util.h
	gcc -c util.c
ssfctp.o: ssfctp.c ssfctp.h
	gcc -c ssfctp.c
clean:
	rm *.o
run:
	./main.exe
Klose:
	(cd "$(shell pwd)/ssfctp_Klose" && gcc SSfctp_main.c ssfctp.c rtime.c readini.c -o main.exe)
runKlose:
	(cd "$(shell pwd)/ssfctp_Klose" && ./main.exe "$(shell pwd)/output_Klose/outfile0.fctp")
NINE := 9
NINETYNINE := 99
FOURSEVENTYNINE := 479
NUMBERS10 := $(shell seq 0 ${NINE})
NUMBERS100 := $(shell seq 0 ${NINETYNINE})
NUMBERS480 := $(shell seq 0 ${FOURSEVENTYNINE})
runKlose10:
	$(foreach var,$(NUMBERS10),(cd "$(shell pwd)/ssfctp_Klose" && ./main.exe "$(shell pwd)/output_Klose/outfile$(var).fctp");)
runKlose100:
	$(foreach var,$(NUMBERS100),(cd "$(shell pwd)/ssfctp_Klose" && ./main.exe "$(shell pwd)/output_Klose/outfile$(var).fctp");)
runKlose480:
	$(foreach var,$(NUMBERS480),(cd "$(shell pwd)/ssfctp_Klose" && ./main.exe "$(shell pwd)/output_Klose/outfile$(var).fctp");)