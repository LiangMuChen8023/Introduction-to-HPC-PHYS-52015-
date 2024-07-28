# Module: Introduction to HPC (PHYS 52015)
**Term:** Summer 2024

**Coursework:** High-performance computing

**Lecturer:** Dr. Christopher Marcotte

---

## Submission
Please submit two PDF files (`part1.pdf` & `part2.pdf`) and two code files (`part1.c` & `part2.c`) in a zip archive.

**Deadlines:** Consult the MISCADA learning and teaching handbook for submission deadlines.

**Plagiarism and collusion:** Students suspected of plagiarism, either of published work or work from unpublished sources, including the work of other students, or of collusion will be dealt with according to Computer Science and University guidelines.

---

## Description
Throughout the course we have considered simple programming problems largely distinct from the typical day-to-day practice of scientific computing. In this assignment, you will experience a sliver of that practice by inheriting a codebase it is your responsibility to make faster while maintaining correctness.

Consider the logistic growth equation,

\[ \frac{du}{dt} = u(1− u), \quad u(t = 0) = u_0, \]

which models the self-limiting growth of a biological population. Under some assumptions, you may represent the state of this system as a map through discrete time (e.g., an annual value) to yield a discrete map or difference equation,

\[ u \leftarrow ru(1− u), \]

where our assumptions of continuity no longer hold. This is an archetypal model of chaotic dynamics for some values of \(r\), and known as the Logistic Map. In the code supplied to you, we step a related difference equation,

\[ u \leftarrow ru(1− u), \]

where \(u\) denotes the local mean of \(u\) in patch size \(\ell^3\) (for \(\ell\) odd) and \(0 \leq r \leq 4\) is a fixed, real-valued, scalar parameter. Around index \((i, j, k)\) the local mean can be written as,

\[ u_{i,j,k} = \ell^{-3} \sum_{(i', j', k') \in N_\ell(i, j, k)} u_{i', j', k'}, \]

such that the neighborhood of index position \((i, j, k)\) is written

\[ N_\ell(i, j, k) = \{ (i', j', k') : \| (i, j, k) - (i', j', k') \|_\infty \leq (\ell - 1)/2 \}, \]

where \(\| \cdots \|_\infty\) is the Chebyshev-distance or \(\infty\)-norm or sup-norm: \(\| x \|_\infty = \max(|x_1|, |x_2|, \cdots)\).

For example, the probability that someone is ill depending on their proximity to others who may be ill, set in a (bizarre) cube environment, and mapped over a discrete interval.

In the provided serial code, you will see that \(u\) is computed by averaging a \(3 \times 3 \times 3\) chunk of the \(u\) array. Of course, at the boundaries of the domain, there are missing elements in the local element \(u_{010}\) has fewer neighbors than \(u_{111}\), and element \(u_{000}\) has fewer neighbors than \(u_{001}\) – which is accounted for by only computing the local patch averages over existing neighbor elements. Your implementation should respect this effect near the physical boundaries, and care should be taken that this is correct for part 2.

### Provided materials
In addition to this document, you are provided with five code files: `serial.c`, `params.h`, `serial.sh`, `part1.sbatch`, and `part2.sbatch`. The file `serial.c` is a serial implementation of the model which you should modify to complete the coursework, and use to verify the correctness of your parallel implementations. The second is a header file, `params.h`, which contains the parameters used in the program. The third file is a simple run script for the serial code. The fourth and fifth files are templates for the compilation and execution of your code on Hamilton; you must ensure that your code functions correctly with these files as provided. Further, example good and bad reports are provided – their content may not be specifically related to this task and any data presented therein is fictionalized. These are meant to serve as guidance and to set expectations for the writing.

### Notes and Warnings
- If you submit new header files with your assignment, it is your responsibility to ensure the submission compiles and runs correctly on Hamilton 8 with the build scripts provided – I will not debug; code which does not run will receive no credit.
- If, for whatever reason, \(u_0 < 0\) or \(u_0 > 1\) then your results will be incorrect; ensure your initial condition is within \(0 \leq u_0 \leq 1\). Likewise, if \(0 \leq u_0 \leq 1\) and \(u\) is ever outside this range then your implementation is incorrect.
- The serial code records the trace of several diagnostic variables, namely the extrema of \(u\), the global mean and variance,

  \[ \mu = N_\infty(i, j, k)^{-1} \sum_{(i, j, k)} u_{i, j, k}, \quad \sigma^2 = N_\infty(i, j, k)^{-1} \sum_{(i, j, k)} |u_{i, j, k} - \mu|^2. \]

  Your code must produce the same outputs as the serial implementation for the same initial conditions.
- You may recognize the local mean computation as a ‘box blur’ or convolution operation. The implementation in the serial code is pedagogical and not optimized. You may be inspired to improve the implementation by refactoring some computations. In your performance investigations, any refactoring done for your parallel codes should also be applied to the serial code and described to ensure you are making a fair performance comparison.
- If you find yourself unable to complete a task, you should address why that is in the report. This may go some way to getting you partial credit, provided you explain your reasoning and demonstrate that you have considered how to do it in some detail.
- Reports which do not demonstrate meeting the learning outcomes or adequately address the content of the reports will receive a low mark; please see the good and bad example reports distributed with this coursework.

### Assessment
In this assessment, you will compile and run a serial three-dimensional logistic map code, and compare its performance against a parallelized version that you will write, and explore this comparison in a report. The serial code is made of six functions: `init`, `dudt`, `step`, `stat`, `write`, and `main`.

#### Part 1: OpenMP
**Code**

The expectations for your parallel implementation are to use OpenMP `#pragmas` to:
- Parallelize the function `init`.
- Parallelize the function `dudt`.
- Parallelize the function `step`.
- Parallelize the function `stat`.

Your code should be in a single C file named `part1.c`. Your code must compile and run with the provided submission script `part1.sbatch` on Hamilton, and produce the same outputs as the serial code for equivalent initial conditions in a file named `part1.dat`.

**Report**

Explain and justify your parallelization strategy, using arguments based on the topics covered in the module. Demonstrate the equivalence of the OpenMP-parallelized program outputs visually – e.g., by comparing the final states, or the statistical outputs in `serial.dat` and `part1.dat`, or something more creative; explain any discrepancies. Investigate the strong scaling of your implementation, reporting scaling results using transferable metrics in your report, and compare to theory-informed expectation. Your report should be no more than one (1) page (plus images), in a file named `part1.pdf`. Additional questions you may wish to consider when writing your report are listed below.

Questions you may wish to consider while writing:
- What options for parallelization are available? Why are some more suitable than others?
- What difficulties arise in the parallelization? Why have you not been asked to parallelize the `write` function?
- Where are the necessary synchronization points? The solution statistics require the generation of a few output numbers from a large array; what patterns are available for this computation?
- How did you avoid data races in your solution? Is your parallelization approach the best option? What alternative approaches could be used?

#### Part 2: MPI
**Code**

Your task is to parallelize the serial code using MPI, by distributing the original problem domain into distinct regions on each process. Particular care should be taken to ensure that your results compile, run, and all outputs are correct (for deterministic inputs). For this part, your code must run correctly with 4 MPI processes. You need not stick to the functions as presented in the serial code. Your implementation should:
- Correctly distribute the allocation of `u` across processes.
- Correctly exchange necessary information across processes.
- Correctly calculate the statistical outputs when `u` is distributed across processes.
- Correctly evaluate the mapping of `u` across processes, taking into account the physical boundaries.

Your code should be a single C file called `part2.c`. Your code should compile and run, producing correct outputs, with the provided submission script (using 4 MPI processes). The outputs should be in a file named `part2.dat`.

**Report**

You should explain and justify your parallelization strategy in a report, demonstrating the equivalence of the MPI-parallelized program outputs visually, and explaining any discrepancies. Your report should also investigate the strong scaling of your implementation and compare it to the OpenMP implementation. Your report should be no more than one (1) page (plus images), in a file named `part2.pdf`.

Questions you may wish to consider while writing:
- How should the domain be decomposed? Why have you chosen your specific decomposition strategy?
- What difficulties arise in the domain decomposition and parallelization?
- How did you ensure that your solution was correct? What test cases were used?
- How does the MPI implementation compare to the OpenMP implementation in terms of performance and scalability?
- What optimizations did you apply to improve performance?
- How did you handle boundary conditions in the distributed computation?

## MyCode
I have successfully implemented the OpenMP and MPI code, and the results are as follows. 
![1722164281977](https://github.com/user-attachments/assets/ffc9ca51-8145-4855-9c91-cef0b43de205)

If you need the code or programming assistance, please contact me on WeChat lqs_8023 or phone 18711839961

