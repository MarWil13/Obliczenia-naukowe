# Numeryczne modelowanie dyfuzji ciepła w pręcie wielomateriałowym

## 1. Cel projektu

Celem projektu jest implementacja i analiza metod numerycznych dla jednowymiarowego równania przewodnictwa ciepła w pręcie złożonym z jednego lub kilku materiałów. Projekt ma charakter obliczeń naukowych: problem fizyczny jest sprowadzony do równania różniczkowego cząstkowego, a następnie rozwiązywany metodami numerycznymi.

Rozważane równanie ma postać:

```text
rho(x) c(x) dT/dt = d/dx( k(x) dT/dx )
```

gdzie:

- `T(x,t)` jest temperaturą,
- `k(x)` jest przewodnością cieplną,
- `rho(x)` jest gęstością,
- `c(x)` jest ciepłem właściwym.

## 2. Założenia modelu

W projekcie przyjęto następujące założenia:

1. Pręt jest jednowymiarowy.
2. Własności materiałowe są stałe w obrębie każdego segmentu.
3. Na końcach pręta zadane są warunki brzegowe Dirichleta.
4. Przewodność na styku dwóch materiałów jest liczona średnią harmoniczną.
5. Siatka przestrzenna jest jednorodna.

## 3. Co może pójść źle, gdy założenia nie są spełnione?

### Zbyt duży krok czasowy w metodzie jawnej

Metoda Eulera jawnego jest stabilna tylko przy kroku czasowym spełniającym warunek CFL. Jeśli `dt` jest zbyt duże, rozwiązanie może zacząć oscylować, przyjmować niefizyczne wartości lub całkowicie się rozbiec.

### Zła średnia na granicy materiałów

Na styku dwóch materiałów nie należy bezrefleksyjnie stosować średniej arytmetycznej przewodności. Dla materiałów o bardzo różnych przewodnościach poprawniejsza jest średnia harmoniczna, bo słabiej przewodzący materiał ogranicza cały przepływ ciepła.

### Zbyt gruba siatka

Jeżeli liczba punktów siatki jest mała, granica między materiałami jest słabo odwzorowana, a profil temperatury może być niedokładny. Błąd przestrzenny maleje po zagęszczeniu siatki.

### Duży krok w metodzie Cranka-Nicolsona

Metoda Cranka-Nicolsona jest stabilna dla równania dyfuzji, ale przy bardzo dużym kroku czasowym i niegładkich danych początkowych może generować niefizyczne oscylacje.

## 4. Użyte algorytmy

### Dyskretyzacja przestrzenna

Dla punktów wewnętrznych siatki zastosowano konserwatywną dyskretyzację różnicową. Strumień ciepła przez granice komórek zależy od przewodności `k_{i+1/2}` i `k_{i-1/2}`.

Operator przestrzenny ma strukturę trójprzekątniową, dlatego jest przechowywany jako macierz rzadka.

### Schemat theta

W czasie zastosowano rodzinę schematów theta:

```text
(I - theta dt L) T^{n+1} = (I + (1-theta) dt L) T^n + wkład warunków brzegowych
```

Szczególne przypadki:

- `theta = 0` — Euler jawny,
- `theta = 1` — Euler niejawny,
- `theta = 0.5` — Crank-Nicolson.

## 5. Złożoność obliczeniowa

Niech `N` oznacza liczbę punktów siatki, a `M` liczbę kroków czasowych.

- Budowa siatki: `O(N)`.
- Budowa operatora: `O(N)`.
- Jeden krok metody jawnej: `O(N)`.
- Metody niejawne: macierz układu jest faktoryzowana raz, a następnie wielokrotnie używana w kolejnych krokach.
- Pamięć: `O(N)` dla pojedynczego profilu temperatury oraz `O(SN)`, gdy zapisujemy `S` profili w historii rozwiązania.

## 6. Alternatywne podejścia

Możliwe alternatywy dla zastosowanego rozwiązania:

1. Metoda elementów skończonych — lepsza dla geometrii bardziej złożonych niż prosty pręt.
2. Adaptacyjny krok czasowy — pozwala zwiększyć wydajność przy zachowaniu dokładności.
3. Niejednorodna siatka przestrzenna — przydatna, gdy chcemy zagęścić punkty przy granicach materiałów.
4. Metody wyższego rzędu w czasie — mogą poprawić dokładność przy gładkich rozwiązaniach.

## 7. Kluczowe fragmenty kodu

### `build_rod`

Funkcja tworzy siatkę przestrzenną i przypisuje punkty do segmentów materiałowych. Dzięki temu można modelować pręt jednorodny albo złożony, np. stal–miedź.

### `harmonic_mean` i `face_conductivity`

Te funkcje odpowiadają za poprawne obliczenie przewodności na granicach między punktami siatki. Jest to szczególnie ważne dla prętów wielomateriałowych.

### `build_operator`

Funkcja buduje rzadki operator przewodzenia ciepła. To centralny element obliczeń, ponieważ opisuje wpływ sąsiednich punktów siatki na zmianę temperatury w czasie.

### `solve_rod_theta`

Jest to uniwersalny solver. Zamiast pisać trzy osobne implementacje, projekt wykorzystuje jedną funkcję z parametrem `theta`.

## 8. Analiza dokładności

Dokładność zależy od dwóch parametrów:

- kroku przestrzennego `dx`,
- kroku czasowego `dt`.

Dla gładkich rozwiązań dyskretyzacja przestrzenna ma typowo drugi rząd dokładności. W czasie Euler jawny i Euler niejawny są metodami pierwszego rzędu, a Crank-Nicolson jest metodą drugiego rzędu.

W testach porównano stan ustalony dla pręta jednorodnego z rozwiązaniem analitycznym, którym jest liniowy profil temperatury.

## 9. Analiza wydajności

Najtańsza w pojedynczym kroku jest metoda jawna, ale wymaga małego `dt`, więc może potrzebować bardzo wielu kroków. Metody niejawne są droższe w jednym kroku, ale pozwalają stosować znacznie większy krok czasowy. Crank-Nicolson stanowi kompromis między stabilnością i dokładnością.

## 10. Wizualizacja wyników

Wizualizacja może obejmować:

- wykres profilu temperatury `T(x,t)` dla kilku chwil czasu,
- porównanie metod dla tego samego problemu,
- wykres błędu względem rozwiązania analitycznego,
- pokazanie niestabilności metody jawnej po przekroczeniu limitu CFL.

W repozytorium można przechowywać notebook prezentacyjny w katalogu `notebooks/`, a wygenerowane wykresy w katalogu `figures/`.

## 11. Praktyczne znaczenie

Model dyfuzji ciepła w pręcie wielomateriałowym jest uproszczonym modelem wielu problemów inżynierskich, np. przewodzenia ciepła przez połączenia metali, izolacje, elementy konstrukcyjne albo warstwy materiałów. Projekt pokazuje, że wybór metody numerycznej, kroku czasowego i sposobu modelowania granic materiałów ma bezpośredni wpływ na stabilność i wiarygodność wyników.
