import logo from './/images/logo.svg';
import './styles/App.css';
import NavBar from './components/NavBar'
import Main from './Main'
import Footer from './components/Footer'

function App() {
  return (
      <div className="app">
          <NavBar />
          <Main />
          <Footer />
      </div>
  );
}

export default App;
