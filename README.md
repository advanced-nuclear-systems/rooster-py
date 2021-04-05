# ROOSTER: Robust Object-Oriented Solver of Transport Equations for Reactor

One file -- one class.

## Input description

The input file should be named `input` and consits of cards (one line -- one card).

`*` : comment.

`solve fuelgrain` : activates calculations of intragranular gas behaviour.

`solve pointkinetics` **REAC_INS** : activates point reactor kinetics calculations using signal `W` as a reactivity input.

`t0` **0.0** : Initial *time*.

