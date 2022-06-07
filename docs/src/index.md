## Panoramica
Una volta creata la **partizione** del piano o dello spazio indotta da un insieme di oggetti solidi descritti dal proprio bordo, quali superfici in 3D o spezzate poligonali in 2D (singoli tratti di retta), si ottengono i domini del piano e dello spazio, indotti proprio da quella partizione.

Trattasi degli **atomi** di un’algebra booleana: questo significa che con qualche calcolo è possibile trasformare ciascuno degli oggetti originali in vettori binari e poter poi svolgere tra di loro tutte le possibili operazioni booleane, complicate o meno che siano (and, or, not, eccetera), e che dunque non avverranno tra oggetti solidi ma tra vettori binari ove l’operazione verrà eseguita termine a termine.

Il lavoro finale del progetto è quello di descrivere, mediante una matrice binaria, come ogni **generatore** è fatto di atomi.

## Strumenti utilizzati

Lo studio e lo sviluppo del codice sono condotti con l’ausilio di alcuni strumenti di seguito elencati:
- **Debugger**, integrato con il plugin di Julia per il text editor *Visual Studio Code*: utilizzato allo scopo di analizzare a tempo di esecuzione, il comportamento del codice con forte riguardo rispetto a come vengono modificati i dati di input durante l’intero flusso;
- **ProfileView**, *libreria esterna*: utilizzata per graficare lo stack delle chiamate ed individuare l’eventuale presenza di colli di bottiglia;
- **BenchmarkTools**, *libreria esterna*: necessaria per eseguire confronti sui tempi di esecuzione delle varie funzioni e dimostrare che le ottimizzazioni applicate hanno portato ad un miglioramento delle prestazioni.

## Obiettivi

Sono stati individuati alcuni punti su cui ci si è soffermati al fine di migliorare, ove possibile, il progetto:
- **Ottimizzazione delle prestazioni** tramite tecniche di parallelizzazione, di calcolo distribuito e altre migliorie con riferimento al libro *Julia High-Performance*;
- Scrittura di una **documentazione** chiara ed esaustiva per migliorare la leggibilità e la comprensibilità delle funzioni principali individuate utilizzando le docstrings del linguaggio Julia;
- Creazione di file di **test** al fine di verificare l’integrità e la corretta interoperabilità del codice con le altri parti di cui il progetto è composto.