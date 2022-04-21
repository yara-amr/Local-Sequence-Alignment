#Merna Hesham Mahmoud
#Haidy Kamal Makram
#Yara Amr Ahmed
#Maram Nasser
#Maram Hatem
#Yasmine Mohamed El-Gazar


from Bio.SubsMat import MatrixInfo
b62 = MatrixInfo.blosum62


def LocalAlignment(Seqence1, Seqence2):

    UserSelect = int(input("Local Alignment For: \n 1. DNA \n 2. Protein\n"))

    if UserSelect == 1:
        DNAseq1Length = len(Seqence1)
        DNAseq2Length = len(Seqence2)

        DNALoopMatrix = [[0 for x in range(DNAseq1Length + 1)] for y in range(DNAseq2Length + 1)]
        for i in range(DNAseq2Length + 1):
            DNALoopMatrix[i][0] = 0
        for j in range(DNAseq1Length + 1):
            DNALoopMatrix[0][j] = 0
        for i in range(1, DNAseq2Length + 1, 1):
            for j in range(1, DNAseq1Length + 1, 1):
                if Seqence1[j - 1] == Seqence2[i - 1]:
                    MatrixScore = 1
                else:
                    MatrixScore = -2
                DNALoopMatrix[i][j] = max(DNALoopMatrix[i - 1][j] - 1, DNALoopMatrix[i][j - 1] - 1, DNALoopMatrix[i - 1][j - 1] + MatrixScore, 0)
        print(DNALoopMatrix[DNAseq2Length][DNAseq1Length])

        Seq1Result = ""
        Seq2Result = ""
        SeqMatch = ""
        i = DNAseq2Length
        j = DNAseq1Length
        while i > 0 and j > 0:
            up = DNALoopMatrix[i - 1][j] - 1
            left = DNALoopMatrix[i][j - 1] - 1
            if Seqence1[j - 1] == Seqence2[i - 1]:
                MatrixScore = 1
            else:
                MatrixScore = -2
            MatrixDiagonal = DNALoopMatrix[i - 1][j - 1] + MatrixScore
            if DNALoopMatrix[i][j] == MatrixDiagonal:
                Seq1Result += Seqence1[j - 1]
                Seq2Result += Seqence2[i - 1]
                if MatrixScore == 1:
                    SeqMatch += "|"
                else:
                    SeqMatch += " "
                i -= 1
                j -= 1
            elif DNALoopMatrix[i][j] == up:
                Seq1Result += "-"
                SeqMatch += " "
                Seq2Result += Seqence2[i - 1]
                i -= 1
            else:
                Seq1Result += Seqence1[j - 1]
                Seq2Result += "-"
                SeqMatch += " "
                j -= 1

        while (i > 0):
            Seq1Result += "-"
            Seq2Result += Seqence2[i - 1]
            SeqMatch += " "
            i -= 1
        while (j > 0):
            Seq1Result += Seqence1[j - 1]
            Seq2Result += "-"
            SeqMatch += " "
            j -= 1

        Seq1Result = Seq1Result[::-1]
        SeqMatch = SeqMatch[::-1]
        Seq2Result = Seq2Result[::-1]

        print(Seq1Result)
        print(SeqMatch)
        print(Seq2Result)

    elif UserSelect == 2:
        Proteinseq1Length = len(Seqence1)
        Proteinseq2Length = len(Seqence2)

        ProteinLoopMatrix = [[0 for x in range(Proteinseq1Length + 1)] for y in range(Proteinseq2Length + 1)]
        for i in range(Proteinseq2Length + 1):
            ProteinLoopMatrix[i][0] = 0
        for j in range(Proteinseq1Length + 1):
            ProteinLoopMatrix[0][j] = 0
        for i in range(1, Proteinseq2Length + 1, 1):
            for j in range(1, Proteinseq1Length + 1, 1):
                if Seqence1[j - 1] == Seqence2[i - 1]:
                    MatrixScore = b62[((Seqence1[j - 1]), (Seqence1[i - 1]))]
                else:
                    MatrixScore = b62[((Seqence1[j - 1]), (Seqence1[i - 1]))]
                ProteinLoopMatrix[i][j] = max(ProteinLoopMatrix[i - 1][j] - 1, ProteinLoopMatrix[i][j - 1] - 1,
                                       ProteinLoopMatrix[i - 1][j - 1] + MatrixScore, 0)
        print(ProteinLoopMatrix[Proteinseq2Length][Proteinseq1Length])

        Seq1Result = ""
        Seq2Result = ""
        SeqMatch = ""
        i = Proteinseq2Length
        j = Proteinseq1Length
        while i > 0 and j > 0:
            up = ProteinLoopMatrix[i - 1][j] - 1
            left = ProteinLoopMatrix[i][j - 1] - 1
            if Seqence1[j - 1] == Seqence2[i - 1]:
                MatrixScore = b62[((Seqence1[j - 1]), (Seqence1[i - 1]))]
            else:
                MatrixScore = b62[((Seqence1[j - 1]), (Seqence1[i - 1]))]
            MatrixDiagonal = ProteinLoopMatrix[i - 1][j - 1] + MatrixScore
            if ProteinLoopMatrix[i][j] == MatrixDiagonal:
                Seq1Result += Seqence1[j - 1]
                Seq2Result += Seqence2[i - 1]
                if MatrixScore == 1:
                    SeqMatch += "|"
                else:
                    SeqMatch += " "
                i -= 1
                j -= 1
            elif ProteinLoopMatrix[i][j] == up:
                Seq1Result += "-"
                SeqMatch += " "
                Seq2Result += Seqence2[i - 1]
                i -= 1
            else:
                Seq1Result += Seqence1[j - 1]
                Seq2Result += "-"
                SeqMatch += " "
                j -= 1

        while (i > 0):
            Seq1Result += "-"
            Seq2Result += Seqence2[i - 1]
            SeqMatch += " "
            i -= 1
        while (j > 0):
            Seq1Result += Seqence1[j - 1]
            Seq2Result += "-"
            SeqMatch += " "
            j -= 1

        Seq1Result = Seq1Result[::-1]
        SeqMatch = SeqMatch[::-1]
        Seq2Result = Seq2Result[::-1]

        print(Seq1Result)
        print(SeqMatch)
        print(Seq2Result)

    else:
        print("Please Select From List")

# calling the function
seq1 = input("please enter the first sequence:")
seq2 = input("please enter the second sequence:")
LocalAlignment(seq1, seq2)