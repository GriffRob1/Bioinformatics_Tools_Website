import './styles/App.css';
import NavBar from './components/NavBar'
import Main from './Main'
import Footer from './components/Footer'
import {useEffect} from "react";
import toolsList from "./data/toolsList";

function App() {

    //initializes favorites for first-time users
    useEffect(() => {
        if (localStorage.getItem('favorites') === null) {
            localStorage.setItem('favorites', '[]')
        }
    }, []);

    return (
        <div className="app">
            <NavBar />
            <Main toolsList={toolsList}/>
            <Footer />
        </div>
    );
}

export default App;
