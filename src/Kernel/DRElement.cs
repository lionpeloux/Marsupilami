using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Marsupilami.Kernel
{
    public enum ElementTypes
    {
        Cable,
        Bar,
        Beam
    }
    public enum CableElementTypes
    {
        Simple,
        Sliding
    }
    public enum BarElementTypes
    {
        Simple,
        CompressionOnly,
        TensionOnly,
    }
    public enum BeamElementTypes
    {
        Simple,
        Torsion
    }

    public abstract class DRElement
    {
        #region ========== CHAMPS ==========
        private int _numNodes;      // NOMBRE DE NOEUDS DE L'ELEMENT
        private int _index;         // INDEX DE l'ELEMENT DANS LE TABLEAU GENERAL
        private ElementTypes _type; // TYPE DE L'ELEMENT
        private string _name;       // NOM OU REFERENCE DE L'ELEMENT

        private double _E;          // MODULE D'YOUNG (Pa)
        private double _S;          // SECTION (m²)
        private double _I;          // INERTIE (m4)
        private double _ES;         // RAIDEUR AXIALE
        private double _EI;         // RAIDEUR EN FLEXION

        private Vector3d[] _Xi;     // VERTEX INITIAUX
        private Vector3d[] _Fext;   // EFFORTS EXTERIEURS
        private double[] _L0;       // LONGUEUR A VIDE

        private double[] _Lt;       // NORME EDGE 1->2
        private Vector3d[] _X12;    // EDGE 1->2
        private Vector3d[] _N12;    // EFFORT NORMAL 1->2
        private Vector3d[] _M12;    // MOMENT 1->2
        private Vector3d[] _F1;     // RESULTANTE DES FORCES AU NOEUD 1 (?)
        private Vector3d[] _F2;     // RESULTANTE DES FORCES AU NOEUD 2 (?)
        

        private double _dt;         // PAS DE TEMPS
        private double _g;          // COEFFICIENT DR
        #endregion

        #region ========== VARIABLES LOCALES ==========
        protected double alpha;     // 0.5*dt*dt
        protected double k;         // RAIDEUR
        protected double lm;        // LUMPED MASS
        #endregion

        #region ========== ACCESSEURS ==========
        public int Count
        {
            get { return _numNodes; }
            protected set { _numNodes = value; }
        }
        public int Index
        {
            get { return _index; }
            protected set { _index = value; }
        }
        public string Name
        {
            get { return _name; }
            protected set { _name = value; }
        }
        public ElementTypes Type
        {
            get { return _type; }
            protected set { _type = value; }
        }

        public double E
        {
            get { return _E; }
            set { _E = value; _ES = _E * _S; _EI = _E * _I; }
        }
        public double S
        {
            get { return _S; }
            set { _S = value; _ES = _E * _S; }
        }
        public double I
        {
            get { return _I; }
            set { _I = value; _EI = _E * _I; }
        }
        public double ES
        {
            get { return _ES; }
            protected set { _ES = value; }
        }
        public double EI
        {
            get { return _EI; }
            protected set { _EI = value; }
        }

        public Vector3d[] Xi
        {
            get { return _Xi; }
            set { _Xi = value; }
        }

        public double[] L0
        {
            get { return _L0; }
            protected set { _L0 = value; }
        }
        public Vector3d[] Fext
        {
            get { return _Fext; }
            protected set { _Fext = value; }
        }

        public double[] Lt
        {
            get { return _Lt; }
            protected set { _Lt = value; }
        }
        public Vector3d[] X12
        {
            get { return _X12; }
            protected set { _X12 = value; }
        }
        public Vector3d[] N12
        {
            get { return _N12; }
            protected set { _N12 = value; }
        }
        public Vector3d[] M12
        {
            get { return _M12; }
            protected set { _M12 = value; }
        }
        public Vector3d[] F1
        {
            get { return _F1; }
            protected set { _F1 = value; }
        }
        public Vector3d[] F2
        {
            get { return _F2; }
            protected set { _F2 = value; }
        }

        public double dt
        {
            get { return _dt; }
            protected set { _dt = value; alpha = 0.5 * dt * dt; }
        }
        public double g
        {
            get { return _g; }
            protected set { _g = value; }
        }
        #endregion

        protected DRElement(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
        {
            // MODELE
            this.Index = index;
            this.Count = Xi.Length;

            // GEOMETRIE INITIALE
            this.Xi = Xi;
            this.L0 = L0;

            // EFFORTS EXTERIEURS
            this.Fext = Fext;

            // RELAX
            this.dt = dt;
            this.g = 1;

            // INSTANCIATION DES TABLEAUX
            Lt = new double[Count - 1];
            X12 = new Vector3d[Count - 1];
            N12 = new Vector3d[Count - 1];
            M12 = new Vector3d[Count];
            F1 = new Vector3d[Count - 1];
            F2 = new Vector3d[Count - 1];
        }
        public abstract override string ToString(); // A UTILIER POUR AFFICHER LES INFOS D'UN ELEMENT

        // METHODES DR ===================================================
        public abstract void Update_R(ref Vector3d[][] Rt, ref Vector3d[][] Fext, ref Vector3d[][] Xt); // Update des efforts intérieurs
        public abstract void Update_LM(ref double[][] LMt);                                            // Update des masses fictives

        // METHODES DATA =================================================
        public Vector3d[] X(ref Vector3d[][] Xt)
        {
            return Xt[Index];
        }
        public Vector3d[] X(Vector3d[][] Xt)
        {
            return Xt[Index];
        }
        public Vector3d[] V(ref Vector3d[][] Vt)
        {
            return Vt[Index];
        }
        public Vector3d[] R(ref Vector3d[][] Rt)
        {
            return Rt[Index];
        }
        public Vector3d[] F(ref Vector3d[][] Fext)
        {
            return Fext[Index];
        }
        public double[] LM(ref double[][] LMt)
        {
            return LMt[Index];
        }

        // METHODES EXTRACTION EFFORTS ===================================
        public abstract double[] N();
        public abstract double[] T();
        public abstract double[] M();
    }

    public abstract class DRCableElement : DRElement
    {
        public CableElementTypes cableElementType;

        private DRCableElement(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
            : base(index, Xi, L0, Fext, dt)
        {
            // MODELE
            this.Type = ElementTypes.Cable;

            // PROPRIETE MECA (DEFAUT)
            this.E = 210e9;
            this.S = 100e-6;
            this.I = 0;
            this.ES = E * S;
            this.EI = 0;
        }
        public override double[] N()
        {
            double[] N = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                N[i] = N12[i].Length;
            }

            return N;
        }
        public override double[] T()
        {
            double[] T = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                T[i] = 0.00;
            }

            return T;
        }
        public override double[] M()
        {
            double[] M = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                M[i] = 0.00;
            }

            return M;
        }

        // STATIC CONSTRUCTORS
        public static DRCableElement CreateSimpleCableElement(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
        {
            DRCableElement ce = new SimpleCable(index, Xi, L0, Fext, dt);
            Rhino.RhinoApp.WriteLine(ce.ToString());
            return ce;
        }
        public static DRCableElement CreateSlindingCableElement(int index, Vector3d[] Xi, double LT0, Vector3d[] Fext, double dt)
        {
            DRCableElement ce = new SlidingCable(index, Xi, LT0, Fext, dt);
            Rhino.RhinoApp.WriteLine(ce.ToString());
            return ce;
        }

        // SUB-CLASS
        private class SimpleCable : DRCableElement
        {
            public SimpleCable(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
                : base(index, Xi, L0, Fext, dt)
            {
                this.cableElementType = CableElementTypes.Simple;
            }
            public override string ToString()
            {
                return "CableElement[SimpleCable] = " + Index;
            }

            public override void Update_R(ref Vector3d[][] Rt, ref Vector3d[][] Fext, ref Vector3d[][] Xt)
            {
                for (int nj = 0; nj < Count; nj++)
                { Rt[Index][nj] = Fext[Index][nj]; }

                for (int nj = 0; nj < Count - 1; nj++)
                {
                    X12[nj] = Xt[Index][nj + 1] - Xt[Index][nj];
                    Lt[nj] = X12[nj].Length;
                }

                // Calcul des efforts dans les éléments et somme des résultantes aux noeuds
                for (int nj = 0; nj < Count - 1; nj++)
                {
                    if (Lt[nj] > L0[nj]) { N12[nj] = ES * (1 / Lt[nj] - 1 / L0[nj]) * X12[nj]; }
                    else { N12[nj] = new Vector3d(0.0, 0.0, 0.0); }

                    Rt[Index][nj] -= N12[nj];
                    Rt[Index][nj + 1] += N12[nj];
                }
            }
            public override void Update_LM(ref double[][] LMt)
            {
                LMt[Index][0] = 0;

                for (int nj = 0; nj < Count - 1; nj++)
                {
                    k = (ES / L0[nj] + g * (N12[nj].Length) / Lt[nj]);
                    lm = alpha * k;
                    LMt[Index][nj] += lm;
                    LMt[Index][nj + 1] = lm;
                }
            }
        }
        private class SlidingCable : DRCableElement
        {
            protected double LT0; // LONGUEUR TOTALE AU REPOS
            protected double LTt; // LONGUEUR TOTALE COURANTE
            protected double Nt;  // EFFORT NORMAL DANS LE CABLE

            public SlidingCable(int index, Vector3d[] Xi, double L0, Vector3d[] Fext, double dt)
                : base(index, Xi, new double[1] { L0 }, Fext, dt)
            {
                this.cableElementType = CableElementTypes.Sliding;

                // ATTENTION
                // le tableau L0 perd son sens pour cet element particulier
                // seule la longueur totale à vide nous intéresse pour trouver l'isotension
                LT0 = L0;
            }
            public override string ToString()
            {
                return "CableElement[SlidingCable] = " + Index;
            }

            public override void Update_R(ref Vector3d[][] Rt, ref Vector3d[][] Fext, ref Vector3d[][] Xt)
            {
                for (int nj = 0; nj < Count; nj++)
                { Rt[Index][nj] = Fext[Index][nj]; }

                LTt = 0;
                for (int i = 0; i < Count - 1; i++)
                {
                    X12[i] = Xt[Index][i + 1] - Xt[Index][i];
                    Lt[i] = X12[i].Length;
                    LTt += Lt[i];
                }

                if (LTt > LT0)
                {
                    Nt = ES * (1 / LTt - 1 / LT0); // EFFORT NORMAL UNIFORME
                    for (int i = 0; i < Count - 1; i++)
                    {
                        N12[i] = Nt * (X12[i] / Lt[i]);
                        Rt[Index][i] -= N12[i];
                        Rt[Index][i + 1] += N12[i];
                    }
                }
                else
                {
                    for (int i = 0; i < Count - 1; i++)
                    {
                        N12[i].X = 0; N12[i].Y = 0; N12[i].Z = 0;
                    }
                }
            }
            public override void Update_LM(ref double[][] LMt)
            {
                LMt[Index][0] = 0;

                for (int i = 0; i < Count - 1; i++)
                {
                    //k = (ES / L0[i] + g * (N12[i].Length) / Lt[i]);
                    k = (ES / (LT0 / (double)Count) + g * (N12[i].Length) / LTt); // RAIDEUR HOMOGENE A VERIFIER
                    lm = alpha * k;
                    LMt[Index][i] += lm;
                    LMt[Index][i + 1] = lm;
                }
            }
        }
    }
    public abstract class DRBarElement : DRElement
    {
        public BarElementTypes barElementType;

        private DRBarElement(int index, Vector3d Xi1, Vector3d Xi2, double L0, Vector3d Fext1, Vector3d Fext2, double dt)
            : base(index, new Vector3d[2] { Xi1, Xi2 }, new double[] { L0 }, new Vector3d[2] { Fext1, Fext2 }, dt)
        {
            // MODELE
            this.Type = ElementTypes.Bar;

            // PROPRIETE MECA (DEFAUT)
            this.E = 25e9;
            this.S = 4e-4;
            this.I = 7.9e-8;
            this.ES = E * S;
            this.EI = E * I;
        }
        public override double[] N()
        {
            double[] N = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                N[i] = N12[i].Length;
            }

            return N;
        }
        public override double[] T()
        {
            double[] T = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                T[i] = (M12[i + 1] - M12[i]).Length / Lt[i];
            }

            return T;
        }
        public override double[] M()
        {
            double[] M = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                M[i] = Math.Max(M12[i].Length, M12[i + 1].Length);
            }

            return M;
        }

        // STATIC CONSTRUCTORS
        public static DRBarElement CreateSimpleBarElement(int index, Vector3d Xi1, Vector3d Xi2, double L0, Vector3d Fext1, Vector3d Fext2, double dt)
        {
            DRBarElement be = new SimpleBar(index, Xi1, Xi2, L0, Fext1, Fext2, dt);
            Rhino.RhinoApp.WriteLine(be.ToString());
            return be;
        }

        // SUB-CLASS
        private class SimpleBar : DRBarElement
        {
            public SimpleBar(int index, Vector3d Xi1, Vector3d Xi2, double L0, Vector3d Fext1, Vector3d Fext2, double dt)
                : base(index, Xi1, Xi2, L0, Fext1, Fext2, dt)
            {
                this.barElementType = BarElementTypes.Simple;
            }
            public override string ToString()
            {
                return "BarElement[SimpleBar] = " + Index;
            }

            public override void Update_R(ref Vector3d[][] Rt, ref Vector3d[][] Fext, ref Vector3d[][] Xt)
            {
                for (int nj = 0; nj < Count; nj++)
                { Rt[Index][nj] = Fext[Index][nj]; }

                for (int i = 0; i < Count - 1; i++)
                {
                    X12[i] = Xt[Index][i + 1] - Xt[Index][i];
                    Lt[i] = X12[i].Length;
                }

                // Calcul des efforts dans les éléments et somme des résultantes aux noeuds
                for (int i = 0; i < Count - 1; i++)
                {
                    N12[i] = ES * (1 / Lt[i] - 1 / L0[i]) * X12[i];
                    Rt[Index][i] -= N12[i];
                    Rt[Index][i + 1] += N12[i];
                }
            }
            public override void Update_LM(ref double[][] LMt)
            {
                LMt[Index][0] = 0;

                for (int i = 0; i < Count - 1; i++)
                {
                    k = (ES / L0[i] + g * (N12[i].Length) / Lt[i]);
                    lm = alpha * k;
                    LMt[Index][i] += lm;
                    LMt[Index][i + 1] = lm;
                }
            }
        }
    }
    public abstract class DRBeamElement : DRElement
    {
        public BeamElementTypes beamElementType;

        private DRBeamElement(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
            : base(index, Xi, L0, Fext, dt)
        {
            // MODELE
            this.Type = ElementTypes.Beam;

            // PROPRIETE MECA (DEFAUT)
            this.E = 25e9;
            this.S = 4e-4;
            this.I = 7.9e-8;
            this.ES = E * S;
            this.EI = E * I;
        }
        public override double[] N()
        {
            double[] N = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                N[i] = N12[i].Length;
            }

            return N;
        }
        public override double[] T()
        {
            double[] T = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                T[i] = (M12[i + 1] - M12[i]).Length / Lt[i];
            }

            return T;
        }
        public override double[] M()
        {
            double[] M = new double[Count - 1];

            for (int i = 0; i < Count - 1; i++)
            {
                M[i] = Math.Max(M12[i].Length, M12[i + 1].Length);
            }

            return M;
        }

        // STATIC CONSTRUCTORS
        public static DRBeamElement CreateSimpleBeamElement(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
        {
            DRBeamElement be = new SimpleBeam(index, Xi, L0, Fext, dt);
            Rhino.RhinoApp.WriteLine(be.ToString());
            return be;
        }

        // SUB-CLASS
        private class SimpleBeam : DRBeamElement
        {
            public SimpleBeam(int index, Vector3d[] Xi, double[] L0, Vector3d[] Fext, double dt)
                : base(index, Xi, L0, Fext, dt)
            {
                this.beamElementType = BeamElementTypes.Simple;
            }
            public override string ToString()
            {
                return "BeamElement[SimpleBeam] = " + Index;
            }

            public override void Update_R(ref Vector3d[][] Rt, ref Vector3d[][] Fext, ref Vector3d[][] Xt)
            {
                for (int nj = 0; nj < Count; nj++)
                { Rt[Index][nj] = Fext[Index][nj]; }

                // CALCUL DES EDGES : VECTEURS ET NORMES
                for (int i = 0; i < Count - 1; i++)
                {
                    X12[i] = Xt[Index][i + 1] - Xt[Index][i];
                    Lt[i] = X12[i].Length;
                }

                // CALCUL DES MOMENTS
                M12[0].X = 0;
                M12[0].Y = 0;
                M12[0].Z = 0;
                for (int i = 1; i < Count - 1; i++)
                {
                    M12[i] = -(2 * EI / (Xt[Index][i + 1] - Xt[Index][i - 1]).Length) * Vector3d.CrossProduct(X12[i], X12[i - 1]) / (Lt[i - 1] * Lt[i]);
                }
                M12[Count - 1].X = 0;
                M12[Count - 1].Y = 0;
                M12[Count - 1].Z = 0;

                // CALCUL DES RESULTANTES
                for (int i = 0; i < Count - 1; i++)
                {
                    N12[i] = ES * (1 / Lt[i] - 1 / L0[i]) * X12[i];
                    F1[i] = Vector3d.CrossProduct(X12[i], M12[i]) / (Lt[i] * Lt[i]);
                    F2[i] = Vector3d.CrossProduct(M12[i + 1], X12[i]) / (Lt[i] * Lt[i]);

                    Rt[Index][i] -= N12[i];
                    Rt[Index][i] -= F1[i];
                    Rt[Index][i] -= F2[i];

                    Rt[Index][i + 1] += N12[i];
                    Rt[Index][i + 1] += F1[i];
                    Rt[Index][i + 1] += F2[i];
                }
            }
            public override void Update_LM(ref double[][] LMt)
            {
                LMt[Index][0] = 0;

                for (int i = 0; i < Count - 1; i++)
                {
                    k = (ES / L0[i] + g * (N12[i].Length + F1[i].Length + F2[i].Length) / Lt[i]);
                    lm = alpha * k;
                    LMt[Index][i] += lm;
                    LMt[Index][i + 1] = lm;
                }
            }
        }
    }
}
