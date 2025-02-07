# bb2-public
BeatBox v. 2

File INSTALL contains the standard compilation and installation
instructions, with plenty of detail suitable for more experienced
programmers.

Briefly: 

- download Beatbox tar ball, say into your $HOME directory,
- unpack the Beatbox tar ball; this will create a new directory, 
- change into that directory and issue the following commands: 

autoreconf -fi
./configure --prefix=$HOME
make
make install

If all works well, this will create BeatBox executables both in the current directory and in your $HOME/bin/ directory; you can move them elsewhere if required.


---
For brief overview, including compilation and installation, please see the file

doc/beatbox_quickstart.html

---
For BeatBox documentation see file

doc/beatbox.html
