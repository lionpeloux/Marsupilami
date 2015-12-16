using Grasshopper.Kernel;
using Marsupilami.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace Marsupilami.GHComponent
{
    public class GHComponent_SimpleElastica : GH_Component
    {
        
        #region DYNAMIC OUTPUTS ENGINE
        // specific code to refresh outputs
        // FLAGS
        public enum LoopFlag
        {
            running,    // process is running
            aborted,    // process has been manually aborted by the user pressing the esc key
            stoped,     // process has been stoped before iteration_max was reached
            ended,      // process has reached iteration_max
        }
        private int iteration_current = 1;
        private int iteration_groupby;      // group of iterations between output refresh
        private int iteration_max;          // total iteration to run
        private bool loop_reset = true;
        private bool loop_reset_cache = true;
        private LoopFlag flag = LoopFlag.ended;
        private GH_Document ghDocument;
        private System.Windows.Forms.Timer trampoline = new System.Windows.Forms.Timer(); // needed to restart solution for intermediate outputs refresh
        // METHODS
        private void documentSolutionEnd(object sender, GH_SolutionEventArgs e)
        {
            ghDocument.SolutionEnd -= documentSolutionEnd;

            // The trampoline gives time to stop the infinitly recursive chain from canvas,
            // as well as mantains the call stack fixed in length
            trampoline.Interval = 1;
            trampoline.Tick += trampolineTick;
            trampoline.Start();
        }
        private void trampolineTick(object sender, EventArgs e)
        {
            trampoline.Tick -= trampolineTick;
            trampoline.Stop();
            this.ExpireSolution(true);
        }
        #endregion

        // FIELDS
        private int index;
        private int N, N_cache;
        private Point3d pt1, pt1_cache, pt2, pt2_cache;
        private double L, L_cache;
        private double E, S, I;
        List<Vector3d> Fext, Fext_cache;
        bool fext_haschanged;
        Vector3d[] Xi;
        double[] Lr;
        

        public DRElement[] elements;
        public List<DRConstraint> constraints;
        public DRRelax solver;
        List<double> Ec;
        public double dt;
        public double Eclim;

        // CONSTRUCTOR
        public GHComponent_SimpleElastica()
            : base("Simple Elastica", "Simple Elastica", "Simple Elastica", "Marsupilami", "Example")
        {
            if (DateTime.Now > MarsupilamiInfo.ExpirationDate) { throw new Exception("Enable to launch the gridshell library. Please contact the author."); }

            Fext = new List<Vector3d>();
            Fext_cache = new List<Vector3d>();
            fext_haschanged = false;

            Ec = new List<double>();
            dt = 1;
            Eclim = 1e-20;
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("{BA202347-F42E-455A-96E0-0110B7BC9BB7}"); }
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }

        // PARAMETERS
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Point 1", "pt1", "First support point", GH_ParamAccess.item);
            pManager.AddPointParameter("Point 2", "pt2", "Second support point", GH_ParamAccess.item);
            pManager.AddNumberParameter("Rest Length", "L0", "Rest length", GH_ParamAccess.item, 10);
            pManager.AddIntegerParameter("Number of Elemnts", "N", "Number of elments in the beam", GH_ParamAccess.item, 21);
            pManager.AddNumberParameter("Young's Modulus", "E", "Young's modulus of the beam's material (Pa)", GH_ParamAccess.item, 25e9);
            pManager.AddNumberParameter("Cross Section Area", "S", "Area of the beam's cross section (m2)", GH_ParamAccess.item, 4.2e-4);
            pManager.AddNumberParameter("Area Moment of Inertia", "I", "Area moment of inertia of the beam's cross section (m4)", GH_ParamAccess.item, 7.9e-8);
            pManager.AddVectorParameter("Fext", "Fext", "external forces", GH_ParamAccess.list);

            pManager.AddIntegerParameter("iteration_max", "N_max", "Total number of iterations to run", GH_ParamAccess.item, 100000);
            pManager.AddIntegerParameter("iterations_groupby", "N_refresh", "Group of iterations between output refresh", GH_ParamAccess.item, 1000);
            pManager.AddBooleanParameter("reset", "reset", "Reset the engine", GH_ParamAccess.item);

            pManager[0].Optional = false;
            pManager[1].Optional = false;
            pManager[2].Optional = false;
            pManager[3].Optional = false;
            pManager[4].Optional = false;
            pManager[5].Optional = false;
            pManager[6].Optional = false;
            pManager[7].Optional = true;
            pManager[8].Optional = false;
            pManager[9].Optional = false;
            pManager[10].Optional = false;
        }
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "");
            pManager.Register_PointParam("Xi", "Xi", "");
            pManager.Register_VectorParam("Xt", "Xt", "");
            pManager.Register_VectorParam("Vt", "Vt", "");
            pManager.Register_DoubleParam("LMt", "LMt", "");
            pManager.Register_VectorParam("Rt", "Rt", "");
            pManager.Register_VectorParam("Fext", "Fext", "");
            pManager.Register_DoubleParam("Ect", "Ect", "");
            pManager.Register_DoubleParam("Nt", "Nt", "");
            pManager.Register_DoubleParam("Tt", "Tt", "");
            pManager.Register_DoubleParam("Mt", "Mt", "");
            pManager.Register_IntegerParam("index", "index", "");
        }

        // SOLVER
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ghDocument = this.OnPingDocument();

            Fext.Clear();

            if (!DA.GetData(0, ref pt1) || !DA.GetData(1, ref pt2)) { return; }
            if (!DA.GetData(2, ref L)) { return; }
            if (!DA.GetData(3, ref N)) { return; }
            if (!DA.GetData(4, ref E)) { return; }
            if (!DA.GetData(5, ref S)) { return; }
            if (!DA.GetData(6, ref I)) { return; }
            if (!DA.GetDataList(7, Fext)) { for (int i = 0; i < N + 1; i++) { Fext.Add(new Vector3d(0, 0, 0)); } }

            if (!DA.GetData(8, ref iteration_max)) { return; }
            if (!DA.GetData(9, ref iteration_groupby)) { return; }
            if (!DA.GetData(10, ref loop_reset)) { return; }

            if (Fext_cache.Count == Fext.Count)
            {
                for (int i = 0; i < Fext_cache.Count; i++)
                {
                    if (Fext[i] != Fext_cache[i])
                    {
                        fext_haschanged = true;
                        break;
                    }
                }
            }

            if (solver == null || loop_reset != loop_reset_cache) // Premier Calcul
            {
                // INIT CACHE TRACKING
                pt1_cache = new Point3d(pt1);
                pt2_cache = new Point3d(pt2);
                L_cache = L;
                N_cache = N;
                Fext_cache.Clear(); fext_haschanged = false;
                for (int i = 0; i < Fext.Count; i++) { Fext_cache.Add(new Vector3d(Fext[i])); }

                Xi = new Vector3d[N + 1];
                for (int i = 0; i < N + 1; i++)
                {
                    Xi[i].X = pt1.X + (pt2.X - pt1.X) * ((double)i / N);
                    Xi[i].Y = pt1.Y + (pt2.Y - pt1.Y) * ((double)i / N);
                    Xi[i].Z = pt1.Z + (pt2.Z - pt1.Z) * ((double)i / N);
                }

                // INIT
                loop_reset_cache = loop_reset;
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();
            }
            else if (pt1.DistanceTo(pt1_cache) != 0 || pt2.DistanceTo(pt2_cache) != 0)
            {
                Xi = new Vector3d[N + 1];
                for (int i = 0; i < N + 1; i++) { Xi[i] = new Vector3d(solver.Xt[0][i]); }
                Xi[0] = new Vector3d(pt1);
                Xi[N] = new Vector3d(pt2);

                // INIT
                pt1_cache = new Point3d(pt1);
                pt2_cache = new Point3d(pt2);
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();
            }
            else if (L != L_cache)
	        {
                Xi = new Vector3d[N + 1];
                for (int i = 0; i < N + 1; i++) { Xi[i] = new Vector3d(solver.Xt[0][i]); }

                // INIT
                L_cache = L;
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();	 
	        }
            else if (N != N_cache)
            {
                Point3d[] pts = new Point3d[N_cache + 1];
                for (int i = 0; i < N_cache + 1; i++) { pts[i] = new Point3d(solver.Xt[0][i]); }

                Curve crv = Curve.CreateInterpolatedCurve(pts,3);                
                crv.DivideByCount(N, true, out pts);

                Xi = new Vector3d[N + 1];
                for (int i = 0; i < N + 1; i++) { Xi[i] = new Vector3d(pts[i]); }

                // INIT
                N_cache = N;
                Fext_cache.Clear(); fext_haschanged = false;
                for (int i = 0; i < Fext.Count; i++) { Fext_cache.Add(new Vector3d(Fext[i])); }
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();	
            }
            else if (fext_haschanged)
            {
                Xi = new Vector3d[N + 1];
                for (int i = 0; i < N + 1; i++) { Xi[i] = new Vector3d(solver.Xt[0][i]); }

                // INIT
                Fext_cache.Clear(); fext_haschanged = false;
                for (int i = 0; i < Fext.Count; i++) { Fext_cache.Add(new Vector3d(Fext[i])); }
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();	        
            }
            else
            {
                if (flag == LoopFlag.running)
                {
                    ghDocument.SolutionEnd += documentSolutionEnd;

                    // GROUP OF ITERATIONS
                    while (flag == LoopFlag.running)
                    {
                        // SINGLE ITERATION
                        iteration_current++;
                        Run();

                        // LOOPSTATE : STOP
                        if (flag == LoopFlag.stoped)
                        {
                            ghDocument.SolutionEnd -= documentSolutionEnd;
                            break;
                        }

                        // LOOPSTATE : END
                        if (iteration_current == iteration_max)
                        {
                            ghDocument.SolutionEnd -= documentSolutionEnd;
                            flag = LoopFlag.ended;
                            break;
                        }

                        // LOOP STATE : ABORT
                        if (GH_Document.IsEscapeKeyDown())
                        {
                            ghDocument.SolutionEnd -= documentSolutionEnd;
                            flag = LoopFlag.aborted;
                            break;
                        }

                        // REFRESH DETECTION
                        if (iteration_current % iteration_groupby == 0) { OnRefresh(); break; }
                    }
                }
            }

            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, flag.ToString());

            DA.SetDataList(0, solver.Info());
            DA.SetDataTree(1, solver.GHTree_Xi());
            DA.SetDataTree(2, solver.GHTree_Xt());
            DA.SetDataTree(3, solver.GHTree_Vt());
            DA.SetDataTree(4, solver.GHTree_LMt());
            DA.SetDataTree(5, solver.GHTree_Rt());
            DA.SetDataTree(6, solver.GHTree_Fext());
            DA.SetDataList(7, Ec);
            DA.SetDataTree(8, solver.GHTree_N());
            DA.SetDataTree(9, solver.GHTree_T());
            DA.SetDataTree(10, solver.GHTree_M());

            DA.SetData(11, index);
            index++;
        }

        private void Init()
        {
            // reset flags, do first iteration, refresh outputs and restart iterative process
            index = 0;
            iteration_current = 1;
            flag = LoopFlag.running;

            Lr = new double[N];
            for (int i = 0; i < N; i++) { Lr[i] = L / N; }

            // ELEMENT
            elements = new DRElement[1] { DRBeamElement.CreateSimpleBeamElement(0, Xi.ToArray(), Lr, Fext.ToArray(), dt) };
            elements[0].E = E;
            elements[0].S = S;
            elements[0].I = I;

            // CONSTRAINTS
            constraints = new List<DRConstraint>();
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(0, 0));
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(0, N));

            // SOLVER
            solver = new DRRelax(elements, constraints);
            solver.Eclim = Eclim;
            Ec.Clear();
            Ec.Add(solver.Ect);
        }
        private void Run()
        {
            solver.Run();
            if (solver.Ect < solver.Eclim) { flag = LoopFlag.stoped; }
            Ec.Add(solver.Ect);
        }
        private void OnRefresh()
        {
           
        }
    }
}
