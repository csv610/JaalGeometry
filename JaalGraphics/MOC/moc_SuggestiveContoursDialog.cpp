/****************************************************************************
** Meta object code from reading C++ file 'SuggestiveContoursDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "SuggestiveContoursDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'SuggestiveContoursDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JSuggestiveContoursDialog_t {
    QByteArrayData data[8];
    char stringdata0[118];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JSuggestiveContoursDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JSuggestiveContoursDialog_t qt_meta_stringdata_JSuggestiveContoursDialog = {
    {
QT_MOC_LITERAL(0, 0, 25), // "JSuggestiveContoursDialog"
QT_MOC_LITERAL(1, 26, 11), // "closeDialog"
QT_MOC_LITERAL(2, 38, 0), // ""
QT_MOC_LITERAL(3, 39, 8), // "setParam"
QT_MOC_LITERAL(4, 48, 13), // "smoothNormals"
QT_MOC_LITERAL(5, 62, 15), // "smoothCurvature"
QT_MOC_LITERAL(6, 78, 20), // "smoothCurvatureDeriv"
QT_MOC_LITERAL(7, 99, 18) // "subdivisionSurface"

    },
    "JSuggestiveContoursDialog\0closeDialog\0"
    "\0setParam\0smoothNormals\0smoothCurvature\0"
    "smoothCurvatureDeriv\0subdivisionSurface"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JSuggestiveContoursDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    0,   46,    2, 0x08 /* Private */,
       5,    0,   47,    2, 0x08 /* Private */,
       6,    0,   48,    2, 0x08 /* Private */,
       7,    0,   49,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JSuggestiveContoursDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JSuggestiveContoursDialog *_t = static_cast<JSuggestiveContoursDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->closeDialog(); break;
        case 1: _t->setParam(); break;
        case 2: _t->smoothNormals(); break;
        case 3: _t->smoothCurvature(); break;
        case 4: _t->smoothCurvatureDeriv(); break;
        case 5: _t->subdivisionSurface(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JSuggestiveContoursDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JSuggestiveContoursDialog.data,
      qt_meta_data_JSuggestiveContoursDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JSuggestiveContoursDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JSuggestiveContoursDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JSuggestiveContoursDialog.stringdata0))
        return static_cast<void*>(const_cast< JSuggestiveContoursDialog*>(this));
    if (!strcmp(_clname, "Ui::SuggestiveContoursDialog"))
        return static_cast< Ui::SuggestiveContoursDialog*>(const_cast< JSuggestiveContoursDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JSuggestiveContoursDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
