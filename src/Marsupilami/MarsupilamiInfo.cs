using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Marsupilami
{
    public class MarsupilamiInfo : GH_AssemblyInfo
    {
        //Fields
        private static DateTime _date_expiration = new DateTime(2016, 10, 30);

        //Properties
        public static DateTime ExpirationDate
        {
            get { return _date_expiration; }
        }

        //Propreties Overrided
        public override string Name
        {
            get { return "Marsupilami"; }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get { return new Guid("{D58B578A-0848-4F1B-A485-1AC61B24BFFD}"); }
        }
        public override string AuthorName
        {
            get { return "Lionel du Peloux"; }
        }
        public override string AuthorContact
        {
            get { return "lionel.dupeloux@gmail.com"; }
        }
    }
}
