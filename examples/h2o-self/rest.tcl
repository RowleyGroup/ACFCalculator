colvarsTrajFrequency 10000

colvar {
    name solute1
    distanceZ {
        main {
            atomnumbers { 1 2 3 }
        }

        ref {
            dummyAtom ( 0.000, 0.000, 0.000 )
        }

        axis (0.0, 0.0, 1.0)
    }
}

harmonic {
    name soluterest1
    colvars  solute1
    centers 0.0
    forceConstant { 10 }
}
