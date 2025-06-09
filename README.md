Metody Numeryczne – Projekt
Projekt realizowany w ramach kursu Metody Numeryczne. Zawiera implementacje popularnych algorytmów numerycznych oraz przykładowe użycia i testy.

Struktura katalogów
MetodyNumeryczne/

├── src/ – Pliki źródłowe biblioteki numerycznej (.cpp)

├── headers/ – Pliki nagłówkowe (.h/.hpp)

├── test/ – Pliki testowe

│ └── tests.cpp

├── examples/ – Przykładowe użycia biblioteki

│ └── Example_Linear.cpp

├── Makefile – Plik make do budowania projektu

└── README.md – Niniejsza dokumentacja


Kompilacja i uruchamianie
Wymagania:

Kompilator zgodny z C++17 (np. g++)

System zgodny z make (Linux/MacOS lub Windows z Git Bash, WSL, MinGW)

Komendy:

make
– Buduje cały projekt (tests + example_linear)

make tests
– Buduje tylko plik testowy

make example_linear
– Buduje tylko przykład Example_Linear

make clean
– Usuwa pliki pośrednie i pliki wykonywalne

Komponenty biblioteki
Biblioteka zawiera implementacje metod:

integral.cpp – metody całkowania numerycznego (np. trapezów, Simpsona)

linear.cpp – rozwiązywanie równań liniowych (np. eliminacja Gaussa)

nonlinear.cpp – rozwiązywanie równań nieliniowych (np. Newton-Raphson)

interp.cpp – interpolacja (np. Lagrange’a, Newtona)

ode.cpp – rozwiązywanie równań różniczkowych zwyczajnych (np. Runge-Kutta)

numlib.cpp – inne funkcje pomocnicze

Testy
Testy znajdują się w pliku:

test/tests.cpp

Uruchamianie:

make tests
./tests

Przykład użycia
Przykładowy program:
examples/Example_Linear.cpp – demonstruje rozwiązywanie układu równań liniowych.

Kompilacja i uruchomienie:

make example_linear
./example_linear

Licencja
Ten projekt służy celom edukacyjnym w ramach kursu akademickiego.
Nie jest przeznaczony do produkcyjnego użytku bez odpowiedniej walidacji.

Autorzy
Krzysztof Mazur, Gabriel Nowak
WiMIP / ITE
Prowadzący: [Siwek Aleksander]
