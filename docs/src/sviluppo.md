In questa sezione vengono presentate tutte le modifiche effettuate al codice, al fine di migliorare il progetto, elencadole nelle diverse fasi di studio.

## Sviluppo studio preliminare
Durante lo studio preliminare, il codice è stato modificato rimuovendo tutte le istruzioni ridondanti o inutili e sono state sistemate molte delle instabilità di tipo emerse utilizzando la macro `@code_warntype`.
Sono stati corretti, inoltre, tutti gli errori riscontrati che non permettevano una corretta esecuzione del codice. Di seguito sono elencate tutte le modifiche effettuate:
- La funzione `cat()`, utilizzata ampiamente nel progetto per concatenare due o più array, presenta ora il parametro `dims`;
- Corretti i tipi delle variabili locali `FEs`, `EVs` e `FVs` nelle funzioni `chainbasis2polygons()` e `chainbasis2solids()`;
- Inizializzata la variabile di ambiente `LC_NUMERIC` con valore `C` come suggerito dalla documentazione del package **Triangulate**.

## Sviluppo studio definitivo
Al fine di introdurre migliorie prestazionali, è stata individuata la possibilità di ricorrere all’utilizzo della parallelizzazione multi-thread, applicata con criterio, evitando di limitarsi all’esclusiva applicazione delle notazioni offerte da Julia (si vedano le annotazioni `@threads`, `@inbounds` e `@inline`).

Si è notato infatti che queste non solo andavano ad offrire risultati inattesi ma alle volte costituivano un vero e proprio blocco rispetto a qualsiasi tentativo di esecuzione del progetto, obbligandoci ad annullare queste modifiche e trovare soluzioni differenti. 

Per garantire che i vari thread possano avere accesso alle stesse informazioni condivise, ma che queste siano accessibili solo da un thread alla volta (**thread safety**) sono state modificate tutte le porzioni di codice non sicure.
A tal riguardo è stata rimossa la funzione `push!()`, dove possibile, dichiarando in anticipo la dimensione degli array o appongiandosi ad array ausiliari.

Sono state parallelizzate funzioni, quali:
- `internalPoints2d()` di `bool2d`;
- `getInternalPoint()` di `bool3d`;
- `chainBasis2Solids()` di `bool3d`;

E’ stata realizzata inoltre una versione specifica della funzione `pointInPolygonClassification()` che serve ad identificare se dati dei punti, questi siano interni o esterni ad un poligono.

Dalla scelta dei punti nelle funzioni `settestpoints()` e `settestpoints2d()`, si è riusciti a semplificare la quantità di casistiche *if-else* da un numero di 15 ad un numero di 4 in una prima iterazione e poi da 4 a 2 in una seconda iterazione. Si è osservato, infatti, che i codici di spigolo ( `c_edge` ) effettivamente determinanti nell'algoritmo sono esclusivamente 3 e 15.

## Test

