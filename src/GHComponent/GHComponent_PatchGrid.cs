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
    public class GHComponent_PatchGrid : GH_Component
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
        private Curve crv_H, crv_V;
        
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
        public GHComponent_PatchGrid()
            : base("Patch Grid", "Patch Grid", "Dynamic relaxation of a patch grid", "Marsupilami", "Example")
        {
            if (DateTime.Now > MarsupilamiInfo.ExpirationDate) { throw new Exception("Enable to launch the gridshell library. Please contact the author."); }

            Ec = new List<double>();
            dt = 1;
            Eclim = 1e-15;
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("{18E628D0-31B8-4F99-8239-4856059CFC0E}"); }
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }

        // PARAMETERS
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("H", "H", "", GH_ParamAccess.item);
            pManager.AddCurveParameter("V", "V", "", GH_ParamAccess.item);

            pManager.AddIntegerParameter("iteration_max", "N_max", "Total number of iterations to run", GH_ParamAccess.item, 100000);
            pManager.AddIntegerParameter("iterations_groupby", "N_refresh", "Group of iterations between output refresh", GH_ParamAccess.item, 1000);
            pManager.AddBooleanParameter("reset", "reset", "Reset the engine", GH_ParamAccess.item);
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
        }

        // SOLVER
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ghDocument = this.OnPingDocument();

            if (!DA.GetData(0, ref crv_H)) { return; }
            if (!DA.GetData(1, ref crv_V)) { return; }

            if (!DA.GetData(2, ref iteration_max)) { return; }
            if (!DA.GetData(3, ref iteration_groupby)) { return; }
            if (!DA.GetData(4, ref loop_reset)) { return; }

            if (solver == null || loop_reset != loop_reset_cache) // Premier Calcul
            {
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
            index++;
        }

        private void Init()
        {
            E = 25e9;
            S = 4.2e-4;
            I = 7.8e-9;

            // reset flags, do first iteration, refresh outputs and restart iterative process
            index = 0;
            iteration_current = 1;
            flag = LoopFlag.running;

            Polyline polyline_H, polyline_V;

            crv_H.TryGetPolyline(out polyline_H);
            crv_V.TryGetPolyline(out polyline_V);

            Lr = new double[N];
            for (int i = 0; i < N; i++) { Lr[i] = L / N; }

            // ELEMENT
            elements = new DRElement[2];

            // Polyline H
            int num = polyline_H.Count;

            Vector3d[] Xi = new Vector3d[num];
            for (int i = 0; i < num; i++)
            {
                Xi[i] = new Vector3d(polyline_H[i]);
            }

            double[] L0 = new double[num - 1];
            for (int i = 1; i < num; i++)
            {
                L0[i-1] = polyline_H[i].DistanceTo(polyline_H[i - 1]);
            }

            Vector3d[] Fext = new Vector3d[num];
            for (int i = 0; i < num; i++)
            {
                Fext[i] = new Vector3d(0, 0, 0);
            }

            elements[0] = DRBeamElement.CreateSimpleBeamElement(0, Xi, L0, Fext.ToArray(), dt);
            elements[0].E = E;
            elements[0].S = S;
            elements[0].I = I;

            // Polyline V
            num = polyline_V.Count;

            Xi = new Vector3d[num];
            for (int i = 0; i < num; i++)
            {
                Xi[i] = new Vector3d(polyline_V[i]);
            }

            L0 = new double[num - 1];
            for (int i = 1; i < num; i++)
            {
                L0[i-1] = polyline_V[i].DistanceTo(polyline_V[i - 1]);
            }

            Fext = new Vector3d[num];
            for (int i = 0; i < num; i++)
            {
                Fext[i] = new Vector3d(0, 0, 0);
            }

            elements[1] = DRBeamElement.CreateSimpleBeamElement(1, Xi, L0, Fext.ToArray(), dt);
            elements[1].E = E;
            elements[1].S = S;
            elements[1].I = I;

            // CONSTRAINTS (APUIS)
            constraints = new List<DRConstraint>();
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(0, 0));
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(0, polyline_H.Count-1));
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(1, 0));
            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(1, polyline_V.Count - 1));


            // CONSTRAINTS (LINK)
            for (int nh = 0; nh < polyline_H.Count; nh++)
            {
                for (int nv = 0; nv < polyline_V.Count; nv++)
                {
                    if (polyline_H[nh].DistanceTo(polyline_V[nv])<0.01)
                    {
                        constraints.Add(DRLinkConstraint.CreateSingleJointConstraint(0, nh, 1, nv));
                        break;
                    }
                }
            }

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
