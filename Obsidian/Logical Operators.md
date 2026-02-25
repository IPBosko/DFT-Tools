A) Unary operator takes one proposition and creates a new one. There are four of them: (i) `NOT` ($\neg$); (ii) `Identity` (ID); (iii) `Tautology` (T); (iv) `Contradiction` ($\perp$).

| $p$ | $\neg p$ | ID $p$ | T $p$ | $\perp p$ |
| --- | -------- | ------ | ----- | --------- |
| $t$ | $f$      | $t$    | $t$   | $f$       |
| $f$ | $t$      | $f$    | $t$   | $f$       |
B) Binary operators. There are $2^4 = 16$ binary operators. Some examples are: `AND` ($\land$), `OR` ($\lor$), `XOR` ($\oplus$), `Implication` ($\implies$), `Equivalence` ($\iff$), etc. 

| $p$ | $q$ | $p\land q$ | $p\lor q$ | $p\oplus q$ | $p\implies q$ | $p\iff q$ |
| --- | --- | ---------- | --------- | ----------- | ------------- | --------- |
| $t$ | $t$ | $t$        | $t$       | $f$         | $t$           | $t$       |
| $t$ | $f$ | $f$        | $t$       | $t$         | $f$           | $f$       |
| $f$ | $t$ | $f$        | $t$       | $t$         | $t$           | $f$       |
| $f$ | $f$ | $f$        | $f$       | $f$         | $t$           | $t$       |

Note: $f\implies t$ is $t$ called "*ex falso quod libet*" -- from a false assumption you can derive whatever you want. 
**Theorem**: $(p\implies q) \iff ((\neg q)\implies(\neg p)).$
**Corollary**: We can prove assertions by way of contradiction. 
**Proof**: Two last columns are equivalent which is the proof. 

| $p$ | $q$ | $\neg p$ | $\neg q$ | $p\implies q$ | $\neg q\implies\neg p$ |
| --- | --- | -------- | -------- | ------------- | ---------------------- |
| $t$ | $t$ | $f$      | $f$      | $t$           | $t$                    |
| $t$ | $f$ | $f$      | $t$      | $f$           | $f$                    |
| $f$ | $t$ | $t$      | $f$      | $t$           | $t$                    |
| $f$ | $f$ | $t$      | $t$      | $t$           | $t$                    |
**Remark I**: Agree on decreasing binding strength in the sequence: $\neg,\land,\lor,\implies,\iff$.
**Remark II**: All higher order operators can be constructed from one single binary operator.