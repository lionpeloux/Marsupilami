using Grasshopper;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace Marsupilami.Kernel
{

    public class DRRelax
    {

        #region ========== CHAMPS ==========
        private double _Eclim = 1e-14;  // ENERGIE CINETIQUE LIMITE
        private double _Ect;            // Ec(t+dt/2)
        private double _Ec0;             // Ec(t-3dt/2)
        private double _Ec1;             // Ec(t-dt/2)

        private double _dt = 1;         // PAS DE TEMPS

        private int numElements;        // NOMBRE TOTAL D'ELEMENTS
        private int numNodes;           // NOMBRE TOTAL DE NOEUDS

        private Vector3d[][] _Xi;       // TABLE DES POSITIONS INITIALES
        private Vector3d[][] _Fext;     // TABLE DES EFFORTS EXTERIEURS

        private Vector3d[][] _Xt;       // TABLE DES POSITIONS
        private Vector3d[][] _Vt;       // TABLE DES VITESSES
        private Vector3d[][] _Rt;       // TABLE DES EFFORTS RESULTANTS
        private double[][] _LMt;        // TABLE DES MASSES FICTIVES [NB : en réalité lmt = dt/lm]       

        private DRElement[] elements;       // TABLEAU DES ELEMENTS
        #endregion

        #region ========== VARIABLES LOCALES ==========
        //
        private double q;                       // FACTEUR D'INTERPOLATION
        // TRACKING
        public int numIteration;                // NOMBRE TOTAL D'ITERATIONS
        public int numPic;                      // NOMBRE DE PIC D'ENERGIE CINETIQUE
        public int numCurrentPicIterations;     // NOMBRE D'ITERATIONS DANS LE PIC EN COURS
        public List<int> iterationHistory;      // HISTORIQUE DES ITERATIONS PAR PIC
        public List<double> chronoHistory;      // HISTORIQUE DES TEMPS PAR PIC
        public Stopwatch chrono;                // TEMPS D'EXECUTION
        // CALCUL VITESSE
        private Vector3d[][] Vtmp;              // TABLE DES VITESSES POUR PIC D'ENERGIE CINETIQUE
        // CONTRAINTES PAR TYPE
        private DRConstraint[] constraints; // TABLEAU DES CONTRAINTES
        private DRKinematicConstraint[] kinematicConstraints;
        private DRMechanicalConstraint[] mechanicalConstraints;
        private DRLinkConstraint[] linkConstraints;
        #endregion

        #region ========== ACCESSEURS ==========
        public int CountElements
        {
            get { return numElements; }
            private set { numElements = value; }
        }
        public int CountNodes
        {
            get { return numNodes; }
            private set { numNodes = value; }
        }
        public double Ect
        {
            get { return _Ect; }
            private set { _Ect = value; }
        }
        public double Ec0
        {
            get { return _Ec0; }
            private set { _Ec0 = value; }
        }
        public double Ec1
        {
            get { return _Ec1; }
            private set { _Ec1 = value; }
        }
        public double Eclim
        {
            get { return _Eclim; }
            set { _Eclim = value; }
        }
        public double dt
        {
            get { return _dt; }
            private set { _dt = value; }
        }
        public Vector3d[][] Fext
        {
            get { return _Fext; }
            private set { _Fext = value; }
        }
        public Vector3d[][] Xi
        {
            get { return _Xi; }
            private set { _Xi = value; }
        }
        public Vector3d[][] Xt
        {
            get { return _Xt; }
            set { _Xt = value; }
        }
        public Vector3d[][] Vt
        {
            get { return _Vt; }
            private set { _Vt = value; }
        }
        public Vector3d[][] Rt
        {
            get { return _Rt; }
            private set { _Rt = value; }
        }
        public double[][] LMt
        {
            get { return _LMt; }
            private set { _LMt = value; }
        }
        public DRElement[] Elements
        {
            get { return elements; }
            private set { elements = value; }
        }
        #endregion

        // CONSTRUCTOR
        public DRRelax(DRElement[] elements, List<DRConstraint> constraints)
        {
            if (DateTime.Now > MarsupilamiInfo.ExpirationDate)
            {
                throw new Exception("Enable to launch the marsupilami library. Please contact the author (lionel.du_peloux@durashell.fr).");
            }
            // TRACKING
            numPic = 0;
            numIteration = 0;
            iterationHistory = new List<int>();
            chronoHistory = new List<double>();
            chrono = new Stopwatch();
            int ni;
            this.elements = elements;
            numElements = elements.Length;
            numNodes = 0;

            DRConstraint.GetConstraints(constraints, ref kinematicConstraints, ref mechanicalConstraints, ref linkConstraints);

            Xi = new Vector3d[numElements][];
            Xt = new Vector3d[numElements][];
            Fext = new Vector3d[numElements][];
            Vt = new Vector3d[numElements][];
            Vtmp = new Vector3d[numElements][];
            Rt = new Vector3d[numElements][];
            LMt = new double[numElements][];

            for (int eli = 0; eli < numElements; eli++)
            {
                ni = elements[eli].Count;
                numNodes += ni;
                Xi[eli] = elements[eli].Xi;
                Fext[eli] = elements[eli].Fext;
                Xt[eli] = (Vector3d[])(elements[eli].Xi.Clone()); // Attention au passage par référence
                Vt[eli] = new Vector3d[ni];
                Vtmp[eli] = new Vector3d[ni];
                Rt[eli] = new Vector3d[ni];
                LMt[eli] = new double[ni];                
            }

            Reset();
        }

        // CORE FUNCTIONS
        private void Init()
        {
            // DO STUFF BEFORE INCREMENTING THE SOLVER
            // - APPLIED DISPLACEMENTS
        }
        private void Reset()
        {
            // RESET PARTICLES VELOCITY TO 0
            chrono.Start();

            // UPDATE NON-CONSTANT EXTERNAL FORCES
            //Update_F();

            // RESET LUMPED MASS & COMPUTE INTERNAL FORCES
            for (int ei = 0; ei < numElements; ei++)
            {
                elements[ei].Update_R(ref _Rt, ref _Fext, ref _Xt);
                elements[ei].Update_LM(ref _LMt);
            }

            // ENFORCE MECHANICAL CONSTRAINTS
            foreach (DRMechanicalConstraint mechanicalConstraint in mechanicalConstraints)
            { mechanicalConstraint.Update_R(ref _Xt, ref _Rt); }

            // TRANSFER NODAL FORCES AND LUMPED MASS BETWEEN LINKED NODES 
            foreach (DRLinkConstraint linkConstraint in linkConstraints)
            { linkConstraint.Transfer_R(ref _Rt); linkConstraint.Transfer_LM(ref _LMt); }

            // RESET VELOCITY (assuming V(0) = 0)
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // V(0 + dt/2) = 1/2 x (dt/lm) x R(0)
                    Vt[ei][nj] = (0.5 * dt / LMt[ei][nj]) * Rt[ei][nj];               
                }
            }

            // ENFORCE VELOCITY CONSTRAINTS
            foreach (DRKinematicConstraint kinematicConstraint in kinematicConstraints)
            { kinematicConstraint.Enforce_V(ref _Xt, ref _Vt); }

            // COMPUTE POSITION
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // Xt(0 + dt) = Xt(0) + dt x V(0 + dt/2)
                    Xt[ei][nj] += dt * Vt[ei][nj];
                }
            }

            // ENFORCE POSITION CONSTRAINTS
            foreach (DRKinematicConstraint kinematicConstraint in kinematicConstraints)
            { kinematicConstraint.Enforce_X(ref _Xt, ref _Vt); }
            foreach (DRLinkConstraint linkConstraint in linkConstraints)
            { linkConstraint.Enforce_X(ref _Xt); }

            // COMPUTE KINETIC ENERGY
            Ect = 0;
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // Ec(0 + dt/2) = S[ 1/2 x lm x V(0 + dt/2)^2 ]
                    Ect += LMt[ei][nj] * Vector3d.Multiply(Vt[ei][nj], Vt[ei][nj]);
                }
            }
            Ect = 0.5 * Ect;    // factorisation de la multiplication par 1/2
            Ec1 = Ect;          // Ec(t -)
            Ec0 = Ect;          // Ec(t -)

            chrono.Stop();
            Rhino.RhinoApp.WriteLine("Ect[" + numPic + "] = " + string.Format("{0:E2}", Ect));
        }
        public void Run()
        {
            // STEP FORWARD
            chrono.Start();

            // UPDATE NON-CONSTANT EXTERNAL FORCES
            //Update_F();

            // COMPUTE INTERNAL FORCES
            for (int ei = 0; ei < numElements; ei++)
            {
                // R(t)
                //elements[ei].Reset_R(ref _Rt, ref _Fext);
                elements[ei].Update_R(ref _Rt, ref _Fext, ref _Xt);
            }

            // ENFORCE MECHANICAL CONSTRAINTS
            foreach (DRMechanicalConstraint mechanicalConstraint in mechanicalConstraints)
            { mechanicalConstraint.Update_R(ref _Xt, ref _Rt); }

            // TRANSFER NODAL FORCES BETWEEN LINKED NODES 
            foreach (DRLinkConstraint linkConstraint in linkConstraints)
            { linkConstraint.Transfer_R(ref _Rt); }

            // COMPUTE VELOCITY
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // V(t + dt/2) = V(t - dt/2) + (dt/lm) x R(t)
                    Vt[ei][nj] += (dt / LMt[ei][nj]) * Rt[ei][nj];
                }
            }

            // ENFORCE VELOCITY CONSTRAINTS
            foreach (DRKinematicConstraint kinematicConstraint in kinematicConstraints)
            { kinematicConstraint.Enforce_V(ref _Xt, ref _Vt); }

            // COMPUTE POSITION
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // X(t + dt) = X(t) + dt x V(t + dt/2)
                    Xt[ei][nj] += dt * Vt[ei][nj];
                }
            }

            // ENFORCE POSITION CONSTRAINTS
            foreach (DRKinematicConstraint kinematicConstraint in kinematicConstraints)
            { kinematicConstraint.Enforce_X(ref _Xt, ref _Vt); }
            foreach (DRLinkConstraint linkConstraint in linkConstraints)
            { linkConstraint.Enforce_X(ref _Xt); }

            // COMPUTE KINETIC ENERGY
            Ect = 0;
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // Ec(t + dt/2) = S[ 1/2 x lm x V(t + dt/2)^2 ]
                    Ect += LMt[ei][nj] * Vector3d.Multiply(Vt[ei][nj], Vt[ei][nj]);
                }
            }
            Ect = 0.5 * Ect;

            if (Ect > Ec1) // KINETIC PIC NOT REACHED
            {
                Ec0 = Ec1; // Ec(t - 3dt/2) <= Ec(t - dt/2)
                Ec1 = Ect; // Ec(t - dt/2) <= Ec(t + dt/2)

                // OLD A VERIFIER entraine (EC0 = EC1)
                // Ec1 = Ect; // Ec(t - dt/2) <= Ec(t + dt/2)
                // Ec0 = Ec1; // Ec(t - 3dt/2) <= Ec(t - dt/2)
            }
            else // KINETIC PIC REACHED
            {
                // X(t*) = f[Ec(t - 3dt/2), Ec(t - dt/2), Ec(t + dt/2)]
                InterpolateEc(Ec0, Ec1, Ect);
                Reset();

                // TRACKING
                numPic++;
                iterationHistory.Add(numIteration);
                chronoHistory.Add(chrono.ElapsedMilliseconds);
            }

            // TRACKING
            numIteration++;
            chrono.Stop();
        }
        private void InterpolateEc(double E0, double E1, double E2)
        {
            // COMPUTE PIC INTERPOLATION
            q = (E1 - E0) / (E0 - 2 * E1 + E2);

            // COMPUTE PREVIOUS VELOCITY
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // Vtmp = V(t-dt/2) = V(t+dt/2) - (dt/lm) x R(t)
                    Vtmp[ei][nj] = Vt[ei][nj] - (dt / LMt[ei][nj]) * Rt[ei][nj];
                }
            }

            // ENFORCE VELOCITY CONSTRAINTS
            foreach (DRKinematicConstraint kinematicConstraint in kinematicConstraints)
            { kinematicConstraint.Enforce_V(ref _Xt, ref Vtmp); }

            // COMPUTE PIC POSITION
            for (int ei = 0; ei < numElements; ei++)
            {
                for (int nj = 0; nj < Xi[ei].Length; nj++)
                {
                    // X(t*) = X(t-dt) + q x dt x V(t-dt/2) = X(t+dt) - dt x V(t+dt/2) + q x dt x V(t-dt/2)
                    Xt[ei][nj] = Xt[ei][nj] - dt * Vt[ei][nj] + q * dt * Vtmp[ei][nj];
                }
            }

            // ENFORCE POSITION CONSTRAINTS
            foreach (DRKinematicConstraint kinematicConstraint in kinematicConstraints)
            { kinematicConstraint.Enforce_X(ref _Xt, ref Vtmp); }
            foreach (DRLinkConstraint linkConstraint in linkConstraints)
            { linkConstraint.Enforce_X(ref _Xt); }

        }

        // SOLVER MONITORING
        public List<string> Info()
        {
            List<string> info = new List<string>();

            info.Add("TOTAL : " + numIteration + " itérations en " + chrono.ElapsedMilliseconds + " ms");
            info.Add(string.Format("COUT : {0:F2} micros/iteration", chrono.ElapsedMilliseconds * 1000 / (double)numIteration));


            if (iterationHistory.Count == 0)
            {
                info.Add("------------------------------------");
            }
            else
            {
                info.Add("NOMBRE DE PICS : " + iterationHistory.Count);
                info.Add("------------------------------------");

                info.Add("P" + (0 + 1) + " : " + iterationHistory[0] + " itérations en " + chronoHistory[0] + " ms");
                for (int i = 0; i < chronoHistory.Count; i++)
                {
                    info.Add("P" + (i + 1) + " : " + iterationHistory[i] + " itérations en " + chronoHistory[i] + " ms");
                }
            }

            return info;
        }

        // IMPLEMENTATION TMP
        private void Update_F()
        {
            //Rhino.RhinoApp.WriteLine("Up");
            int index = 0;
            double p = -15;
            Vector3d t, w, next;
            int n = Xi[0].Length;

            Vector3d[] Ftmp = new Vector3d[n];
            Vector3d[] X12 = new Vector3d[n - 1];

            for (int i = 1; i < n - 1; i++)
            {
                t = Xt[index][i + 1] - Xt[index][i - 1];
                t.Unitize();

                w = Vector3d.CrossProduct(Xt[index][i] - Xt[index][i - 1], Xt[index][i + 1] - Xt[index][i]);
                w.Unitize();

                next = Vector3d.CrossProduct(t, w);
                next.Unitize();

                if (next.Length == 0)
                {
                    Ftmp[i] = -p * Vector3d.ZAxis;
                }
                else
                {
                    Ftmp[i] = -p * next;
                }
            }

            //Rhino.RhinoApp.WriteLine(Ftmp[i].ToString());
            _Fext[index] = Ftmp;
        } // FORCE SUIVEUSES

        // OUTPUTS FOR GRASSHOPPER COMPONENTS
        public DataTree<Vector3d> GHTree_Xt()
        {
            DataTree<Vector3d> tree = new DataTree<Vector3d>();
            int branch = 0;
           

            foreach (DRElement element in elements)
            {
                tree.AddRange(element.X(ref _Xt), new GH_Path(branch));
                branch++;
            }
            return tree;
        }
        public DataTree<Vector3d> GHTree_Xi()
        {
            DataTree<Vector3d> tree = new DataTree<Vector3d>();
            int branch = 0;

            foreach (DRElement element in elements)
            {
                tree.AddRange(element.Xi, new GH_Path(branch));
                branch++;
            }

            return tree;
        }
        public DataTree<Vector3d> GHTree_Vt()
        {
            DataTree<Vector3d> tree = new DataTree<Vector3d>();
            int branch = 0;

            foreach (DRElement element in elements)
            {
                tree.AddRange(element.V(ref _Vt), new GH_Path(branch));
                branch++;
            }
            return tree;
        }
        public DataTree<Vector3d> GHTree_Rt()
        {
            DataTree<Vector3d> tree = new DataTree<Vector3d>();
            int branch = 0;

            foreach (DRElement element in elements)
            {
                tree.AddRange(element.R(ref _Rt), new GH_Path(branch));
                branch++;
            }
            return tree;
        }
        public DataTree<Vector3d> GHTree_Fext()
        {
            DataTree<Vector3d> tree = new DataTree<Vector3d>();
            int branch = 0;

            foreach (DRElement element in elements)
            {
                tree.AddRange(element.F(ref _Fext), new GH_Path(branch));
                branch++;
            }
            return tree;
        }
        public DataTree<GH_Number> GHTree_LMt()
        {
            DataTree<GH_Number> tree = new DataTree<GH_Number>();
            GH_Path path;
            int branch = 0;

            foreach (DRElement element in elements)
            {
                path = new GH_Path(branch);
                foreach (double LMi in element.LM(ref _LMt))
                {
                    tree.Add(new GH_Number(LMi), path);
                }
                branch++;
            }
            return tree;
        }
        public DataTree<GH_Number> GHTree_N()
        {
            DataTree<GH_Number> tree = new DataTree<GH_Number>();
            GH_Path path;
            int branch = 0;

            foreach (DRElement element in elements)
            {
                path = new GH_Path(branch);
                foreach (double Ni in element.N())
                {
                    tree.Add(new GH_Number(Ni), path);
                }
                branch++;
            }
            return tree;
        }
        public DataTree<GH_Number> GHTree_T()
        {
            DataTree<GH_Number> tree = new DataTree<GH_Number>();
            GH_Path path;
            int branch = 0;

            foreach (DRElement element in elements)
            {
                path = new GH_Path(branch);
                foreach (double Ti in element.T())
                {
                    tree.Add(new GH_Number(Ti), path);
                }
                branch++;
            }
            return tree;
        }
        public DataTree<GH_Number> GHTree_M()
        {
            DataTree<GH_Number> tree = new DataTree<GH_Number>();
            GH_Path path;
            int branch = 0;

            foreach (DRElement element in elements)
            {
                path = new GH_Path(branch);

                foreach (double Mi in element.M())
                {
                    tree.Add(new GH_Number(Mi), path);
                }
                branch++;
            }
            return tree;
        }
    }
}
