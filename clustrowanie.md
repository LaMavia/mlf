# uniprot

advanced search, trmd wyszukać, kliknąć w jakiej interpro i z nich wybrać domeny (biolodzy) do advanced search w uniprotie (pod and-ie/or-ze), (może to taxonomy). Patrzeć, czy są od bakterii. Unreviewed mają automatyczne klasyfikacje, więc mogą być złe.

# alignment sekwencyjny - mafft (da się zainstalować)

daje najmniej gapów między alignmentem. Nawet na podstawowych parametrach daje całkiem dobre wyniki. Wyrzuci uliniowione sekwencje.

# clustrowanie sekwencyjne

cdhit (da się zainstalować łatwo, albo zbuildować). intersują -i,-o, sequence id threshold (jak podobne, żeby w tym samym clustrze). zwraca 2 typy plików: zbiór sekwaencji z wybranymi reprezentantami; csltr (jaki podział na klastry: \* - reprezentant, n% - podobieństwo)

# clustrowania strukturalnego

pymol (nie wiemy do końca jaka komenda, ale jak będzie podobieństwo/energie, to wszelkie algorytmy do clustrowania). Da się na colabie używać.

# sequence logo (wizualizacja uliniowienia)

logo/aliview (konwertowanie między formatami alignmentu, ale openbabel też to robi) to tylko wizualizacja: jeśli wykonamy alignment sekwencyjny, to do sprawdzenia sobie.

Jak widać wizualnie: są pozycje, które są lepiej zachowane (conserved). W logo gapy - puste. Można też na frequency ustawić.

Miejsce aktywne w białku musi być zachowane. W przypadku trmd jest to mniej więcej ta sama grupa aminokwasów.

Skąd wiedzieć, które residues to miejsce aktywne: pozycja się nie zmienia, a ta ungapped może. Cechą alignmentu jest to, że sekwencje są równej długości.

Szukamy miejsca aktywnego w alignmencie, liczymy ungapped i wtedy wimy, które selekcje wybrać a pymolu.

Macież blosum.

Wystarczy alignment sekwencyjny i nałożenie go na strukturę w przypadku trmd (chyba™)

Profile sekwencyjne (pewnie niepotrzebne) ~ logo + prawdopodobieństwo danego aminokwasu i niektóre programy tego wymagają (Hidden Markov models, ncbi (paczka), może blast (lokalnie, a nie internety; gdzieś podczas searchu i da się odzyskać))

Podpowiedź potencjalna: ligandMPNN przewidywanie relacji między ligandami i białkami (**ale nowa rzecz, więc nie wiemy, czy będzie ok. lepiej kontynuować pewnie swoje**)
