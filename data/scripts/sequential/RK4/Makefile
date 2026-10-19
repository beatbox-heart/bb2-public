BBXS=bb2s
rk4test_%.txt: rk4test.bbs
	@$(BBXS) $+ $* -mute
ALL=10 30 100 300 1000 3000 10000 30000
ALLTXT=${foreach n,$(ALL),rk4test_$n.txt}
rk4test.txt: $(ALLTXT)
	cat $+ | sort -n > $@
rk4test.eps: rk4test.gnu rk4test.txt
	gnuplot $<
