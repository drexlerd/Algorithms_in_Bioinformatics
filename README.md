# Algorithms_in_Bioinformatics

Interfaces and materials for the programming course "Algorithms in Bioinformatics"

Abstract base classes for the algorithms are located in directory prakt.
Minimal examples to get you started with your implementations are located in the base directory. Examples for pytest unittests are located in directory test. Calling ```pytest``` from the base directory will find and execute these tests. Please keep the script and class names, they are used by the tests to import your derived classes.

The base classes provide a common interface for unittests. You have to implement commandline parameter parsing before calling the run functions and also format the returned results for display on the screen.
