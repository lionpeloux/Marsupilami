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
    public class GHComponent_SimpleCable : GH_Component
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
        private int N;
        private Point3d pt1, pt2;
        private double L;
        private double E, S;
        List<Vector3d> Fext;
        Vector3d[] Xi;
        double[] Lr;
        

        public DRElement[] elements;
        public List<DRConstraint> constraints;
        public DRRelax solver;
        List<double> Ec;
        public double dt;
        public double Eclim;

        // CONSTRUCTOR
        public GHComponent_SimpleCable()
            : base("Simple Cable", "Simple Cable", "Simple Cable", "Marsupilami", "Example")
        {
            if (DateTime.Now > MarsupilamiInfo.ExpirationDate) { throw new Exception("Enable to launch the gridshell library. Please contact the author."); }

            Fext = new List<Vector3d>();

            Ec = new List<double>();
            dt = 1;
            Eclim = 1e-16;
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("{72286637-905E-47DC-B1D2-2D8AA615BA3B}"); }
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
            pManager[6].Optional = true;
            pManager[7].Optional = false;
            pManager[8].Optional = false;
            pManager[9].Optional = false;
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
            if (!DA.GetDataList(6, Fext)) { for (int i = 0; i < N + 1; i++) { Fext.Add(new Vector3d(0, 0, 0)); } }

            if (!DA.GetData(7, ref iteration_max)) { return; }
            if (!DA.GetData(8, ref iteration_groupby)) { return; }
            if (!DA.GetData(9, ref loop_reset)) { return; }


            if (solver == null || loop_reset != loop_reset_cache) // Premier Calcul
            {
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
            elements = new DRElement[1] { DRCableElement.CreateSimpleCableElement(0, Xi.ToArray(), Lr, Fext.ToArray(), dt) };
            elements[0].E = E;
            elements[0].S = S;
            elements[0].I = 0;

            // CONSTRAINTS
            constraints = new List<DRConstraint>();
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(0, 0));
            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(0, N));
            constraints.Add(DRMechanicalConstraint.CreateSpringConstraint(0, N, Xi[N], 10));

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
